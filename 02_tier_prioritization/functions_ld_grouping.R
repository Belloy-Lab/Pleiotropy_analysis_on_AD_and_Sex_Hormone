# =============================================================================
# Helper functions for ld_grouping.R
#
# All paths, the PLINK executable and the genotype reference panel are passed
# in as arguments so that this file contains no machine-specific configuration.
# =============================================================================

# Load a pleioFDR result object for one trait combination.
read_pleiofdr_result <- function(trait, ad_stratum, results_dir) {
  path <- file.path(results_dir,
                    paste0("results_", trait, "_", ad_stratum, "_allexclude"),
                    "result.mat")
  return(R.matlab::readMat(path))
}

# Keep the vectors of length n_ref_snps, bind the LD reference table and drop
# SNPs with missing FDR values.
process_dataframe <- function(data, ref_df, n_ref_snps = 9545380) {
  long_columns <- sapply(data, function(x) length(x) == n_ref_snps)
  data <- as.data.frame(data[long_columns])
  data <- cbind(ref_df, data)
  data <- data[!sapply(data$fdrmat, is.na), ]
  data$fdrmat <- as.numeric(data$fdrmat)
  return(as.data.frame(data))
}

# Read the conjFDR loci table that pleioFDR writes for one trait combination.
read_conjfdr_loci_csv <- function(trait, ad_stratum, results_dir) {
  converted_suffix <- ad_stratum
  if (grepl("ADMeta_female", ad_stratum)) {
    converted_suffix <- "ADMetaF"
  } else if (grepl("ADMeta_male", ad_stratum)) {
    converted_suffix <- "ADMetaM"
  } else if (grepl("ADMeta_all", ad_stratum)) {
    converted_suffix <- "ADMetaA"
  }

  path <- file.path(results_dir,
                    paste0("results_", trait, "_", ad_stratum, "_allexclude"),
                    paste0(converted_suffix, "_", trait, "_conjfdr_0.05_loci.csv"))
  return(data.table::fread(path))
}

# Build the per-locus FDR matrix across all trait combinations.
get_sug_loci_fdr <- function(sug_loci_list, data_list) {
  results_list <- list()

  # Loop through every SNP in sug_loci_list
  for (loci_name in names(sug_loci_list)) {
    cat("Processing loci:", loci_name, "\n")
    loci_data <- sug_loci_list[[loci_name]]

    # Loop through every row in loci_data
    for (i in 1:nrow(loci_data)) {
      temp = paste0(i, "/", nrow(loci_data))
      cat("Processing loci:", loci_name, " ", temp, "\n")
      snp <- loci_data[i,]$SNP
      chr <- loci_data[i,]$CHR
      bp <- loci_data[i,]$BP

      # Find fdrmat value for the SNP in every dataset in data_list
      fdrmat_values <- lapply(data_list, function(x) {
        if (is.null(x)) {
          return(NA)
        }
        idx <- which(x$SNP == snp)
        if (length(idx) == 0) {
          return(NA)
        }
        return(x[idx,]$fdrmat)
      })

      # Combine the results into a single row
      result_row <- c(snp, chr, bp, unlist(fdrmat_values))

      # Add the row to the results list
      results_list <- rbind(results_list, result_row)
    }

    # Print progress
    cat("Processed loci:", loci_name, "\n")
  }

  # Convert results list to data frame
  results_df <- as.data.frame(results_list, stringsAsFactors = FALSE)

  # Add column names for the first three columns
  colnames(results_df) <- c("SNP", "CHR", "BP", paste0(names(data_list)))

  # Convert columns to appropriate data types
  results_df$SNP <- sapply(results_df$SNP, as.character)
  results_df$CHR <- sapply(results_df$CHR, as.numeric)
  results_df$BP <- sapply(results_df$BP, as.numeric)

  # Sort the results by chromosome and position
  results_df <- results_df[order(results_df$CHR, results_df$BP), ]
  results_df <- unique(results_df)

  return(results_df)
}

# Write the SNP list in the BED-like format used by the UCSC liftOver step.
create_bed_file <- function(AD_data, output_file) {

  AD_data$SNP = paste(AD_data$CHR, AD_data$BP, AD_data$A1, AD_data$A2, sep=":")

  fwrite(AD_data, file = paste0(output_file,"_all",".txt"), quote=F, row.names=F, col.names=T, sep='\t')
  AD_data_bed <- AD_data %>%
    dplyr::mutate(start = BP, end = BP) %>%
    dplyr::select(SNP, CHR, start, end)

  fwrite(AD_data_bed, file = paste0(output_file,"_bed",".txt"), sep = "\t", row.names = FALSE, quote = FALSE)
  return(AD_data_bed)
}

# Lift the loci over from hg19 to hg38 with the UCSC liftOver executable.
# `liftover_dir` must contain `liftOver` and `hg19ToHg38.over.chain.gz`.
perform_liftOver <- function(AD_data, liftover_dir, out_fold) {
    chain = "hg19ToHg38.over.chain.gz"

    bed = AD_data

    bed <- as.data.table(bed)
    setnames(bed, c("SNP", "CHR", "start", "end"), c("NAME", "CHR", "START", "END"))
    bed = dplyr::select(bed, CHR, START, END, NAME)

    bed[, CHR := as.character(CHR)]
    bed[, CHR := paste0("chr", CHR)]

    bed[, END := END + 1]

    bed[, START := as.integer(START)]
    bed[, END := as.integer(END)]

    bed_file = paste0(out_fold, "_to_lift_hg19.txt")

    fwrite(bed, bed_file, col.names=F, row.names=F, quote=F, sep='\t')

    liftOver = paste(liftover_dir, 'liftOver', sep="")

    lifted = paste(out_fold, "_lifted_hg38.txt", sep="")
    not_lifted = paste(out_fold, "_not_lifted_hg19tohg38.txt", sep="")

    chain_file = paste(liftover_dir, chain, sep="")

    cmd = paste(liftOver, bed_file, chain_file, lifted, not_lifted, sep=" ")
    system(cmd)

    data_37 = fread(paste0(out_fold,"_all.txt"), sep="\t", header=T)
    bedl = fread(lifted, sep="\t", header=F)

    data_38 = merge(data_37, bedl, by.x = "SNP",by.y = "V4")
    #remove rows where chromosome values do not match
    data_38$V1 = gsub("chr", "", data_38$V1)
    pos = which(data_38$CHR != data_38$V1)
    if (any(pos)) {
        data_38 = data_38[-pos, ]
    }
    data_38 <- data_38 %>%
              dplyr::select(-CHR, -BP) %>%
              dplyr::rename(CHR = V1, BP = V2) %>%
              dplyr::arrange(CHR, BP) %>%
              dplyr::mutate(CHR = gsub("chr", "", CHR)) %>%
              dplyr::select(-V3,-SNP)

    fwrite(data_38, file = paste0(out_fold,"_38.txt"), quote=F, row.names=F, col.names=T, sep='\t')
    return(data_38)
}

# r2 between two SNPs, estimated with PLINK from the genotype reference panel.
# `bfile_prefix` is the PLINK binary prefix without the chromosome number.
get_r2 <- function(snp1, snp2, chr, plink_exe, keep_file, bfile_prefix, work_dir) {
  print(paste("Processing SNP pair:", snp1, snp2, "on chromosome", chr))
  pl <- plink_exe
  in_file_chr <- paste0(bfile_prefix, chr)
  pair_prefix <- file.path(work_dir, "snp_pair")
  pair_txt <- file.path(work_dir, "snp_pair.txt")

  print("Writing SNP pair to temporary file...")
  fwrite(data.frame(snp1, snp2), pair_txt, quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")

  print("Running PLINK command to generate temporary data...")
  command <- paste(pl,
                   "--bfile", in_file_chr,
                   "--allow-no-sex",
                   "--extract", pair_txt,
                   "--keep-allele-order",
                   "--make-bed --out", paste0(pair_prefix, ".data"),
                   sep = " ")
  system(command)

  print("Running PLINK command to calculate r2 matrix...")
  command <- paste(pl,
                   "--bfile", paste0(pair_prefix, ".data"),
                   "--allow-no-sex",
                   "--keep-allele-order",
                   "--keep", keep_file,
                   "--r2 yes-really --matrix",
                   "--out", pair_prefix,
                   sep = " ")
  system(command)

  print("Reading LD matrix...")
  LD <- fread(paste0(pair_prefix, ".ld"))
  LD <- as.matrix(LD)

  temp_files <- c(paste0(pair_prefix, ".ld"),
                  paste0(pair_prefix, ".data.bed"),
                  paste0(pair_prefix, ".data.bim"),
                  paste0(pair_prefix, ".data.fam"),
                  paste0(pair_prefix, ".data.log"),
                  paste0(pair_prefix, ".log"),
                  paste0(pair_prefix, ".data.nosex"),
                  paste0(pair_prefix, ".nosex"),
                  pair_txt)

  if (nrow(LD) < 2) {
    print("LD matrix is empty, cleaning up temporary files...")
    suppressWarnings(file.remove(temp_files))
    return(0)
  }

  r2 <- LD[1, 2]
  print(paste("Calculated r2 value:", r2))
  print("Cleaning up temporary files...")
  suppressWarnings(file.remove(temp_files))

  return(r2)
}

# Group the conjFDR loci into LD groups (one group per independent locus).
get_ld_groups <- function(sug_loci_results, plink_exe, keep_file, bfile_prefix, work_dir,
                          ld_threshold = 0.01) {
  print("Starting LD group calculation...")
  num_snps <- nrow(sug_loci_results)
  print(paste("Number of SNPs:", num_snps))
  ld_matrix <- matrix(0, nrow = num_snps, ncol = num_snps, dimnames = list(sug_loci_results$SNP, sug_loci_results$SNP))

  for (i in 1:(num_snps - 1)) {
    for (j in (i + 1):num_snps) {
      if (sug_loci_results$CHR[i] == sug_loci_results$CHR[j] & abs(sug_loci_results$BP[i] - sug_loci_results$BP[j]) < 1000000) {
        print(paste("Calculating LD for SNPs:", sug_loci_results$SNP[i], sug_loci_results$SNP[j]))
        ld_result <- get_r2(sug_loci_results$SNP[i], sug_loci_results$SNP[j], sug_loci_results$CHR[i],
                            plink_exe = plink_exe,
                            keep_file = keep_file,
                            bfile_prefix = bfile_prefix,
                            work_dir = work_dir)
        ld_matrix[i, j] <- ld_result
        ld_matrix[j, i] <- ld_result
      }
    }
  }

  print("Grouping SNPs based on LD threshold...")
  threshold <- ld_threshold
  groups <- list()
  group_index <- 1
  assigned <- rep(FALSE, num_snps)

  for (i in 1:num_snps) {
    for (j in (i + 1):num_snps) {
      if (j > num_snps) next
      if (ld_matrix[i, j] > threshold) {
        if (!any(unlist(groups) %in% c(i, j))) {
          groups[[group_index]] <- c(i, j)
          assigned[c(i, j)] <- TRUE
          group_index <- group_index + 1
        } else {
          for (k in 1:length(groups)) {
            if (any(groups[[k]] %in% c(i, j))) {
              groups[[k]] <- unique(c(groups[[k]], i, j))
              assigned[c(i, j)] <- TRUE
            }
          }
        }
      }
    }
    if (!assigned[i]) {
      groups[[group_index]] <- c(i)
      assigned[i] <- TRUE
      group_index <- group_index + 1
    }
  }

  print("Assigning groups to SNPs...")
  sug_loci_results$Group <- NA
  for (i in 1:length(groups)) {
    sug_loci_results$Group[groups[[i]]] <- i
  }

  print("Ordering results by group...")
  sug_loci_results <- sug_loci_results[order(sug_loci_results$Group), ]
  sug_loci_results <- unique(sug_loci_results)

  print("LD group calculation complete!")
  return(sug_loci_results)
}
