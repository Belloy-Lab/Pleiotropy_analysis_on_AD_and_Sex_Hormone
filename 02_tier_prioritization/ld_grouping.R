# =============================================================================
# 02_tier_prioritization / ld_grouping.R
#
# Step 1 of the tier-based prioritization:
#   1. collect the conjFDR results of every AD x trait combination
#   2. build the per-locus FDR matrix across all trait combinations
#   3. add allele information from the pleioFDR LD reference
#   4. lift the loci over from hg19 to hg38
#   5. group loci that are in LD (r2 > ld_group_threshold) into one locus group
#   6. write the long-format table of locus/trait pairs with conjFDR < 0.05
#
# Expected layout of the pleioFDR output folder (pleiofdr_results_dir):
#   results_<trait>_<AD stratum>_allexclude/result.mat
#   results_<trait>_<AD stratum>_allexclude/<ADstratumShort>_<trait>_conjfdr_0.05_loci.csv
# with <AD stratum> in {ADMeta_female, ADMeta_male, ADMeta_all} and
# <ADstratumShort> in {ADMetaF, ADMetaM, ADMetaA}.
# =============================================================================

# =========================
# USER CONFIGURATION
# =========================
pleiofdr_results_dir <- "path/to/pleiofdr_output/"
ref_file             <- "path/to/reference/9545380.ref"
liftover_dir         <- "path/to/liftover/"                          # UCSC liftOver executable + hg19ToHg38.over.chain.gz
plink_exe            <- "plink1.9"
ld_keep_file         <- "path/to/reference/unrelated_individuals.txt"
ld_bfile_prefix      <- "path/to/reference_genotypes/EUR_reference_panel_ch"  # chromosome number is appended
ld_work_dir          <- "path/to/work/ld_matrix/"
output_dir           <- "path/to/output/"
n_ref_snps           <- 9545380                       # number of SNPs in the pleioFDR LD reference
ld_group_threshold   <- 0.01                          # r2 threshold used to group the conjFDR loci

# Traits included in the locus table.
#   label : label used in the pleioFDR result folder / file names
#   sex   : which AD stratum the trait is conditioned on ("female" or "male";
#           each trait is additionally conditioned on the non-stratified AD GWAS)
traits <- data.frame(
  label = c("Meno", "Mena", "voice", "BTF", "BTM", "TTF", "TTM",
            "SBFBMI", "SBMBMI", "SBF", "SBM"),
  sex   = c("female", "female", "male", "female", "male", "female", "male",
            "female", "male", "female", "male"),
  stringsAsFactors = FALSE
)
# =========================

library(data.table)
library(dplyr)
library(R.matlab)
library(reshape2)

source("functions_ld_grouping.R")

#-------------------Read Reference Data--------------------
ref_df <- fread(ref_file)

#----------------Initialize Lists to Store Data--------------------
data_list <- list()
sug_loci_list <- list()

#---Loop through datasets and AD combinations to process .mat and .csv files---
for (dataset in traits$label) {

  if (traits$sex[traits$label == dataset] == "female") {
    ad_combinations <- c("ADMeta_female", "ADMeta_all")
  } else {
    ad_combinations <- c("ADMeta_male", "ADMeta_all")
  }

  for (ad_comb in ad_combinations) {

    #----------------Process .mat files--------------------
    variable_name_mat <- paste0(dataset, "_", ad_comb, "_mat")
    raw_data_mat <- read_pleiofdr_result(dataset, ad_comb, pleiofdr_results_dir)
    processed_data_mat <- process_dataframe(raw_data_mat, ref_df, n_ref_snps = n_ref_snps)

    #----------------Sort processed_data_mat by CHR and BP--------------------
    processed_data_mat <- processed_data_mat[order(processed_data_mat$CHR, processed_data_mat$BP), ]

    data_list[[variable_name_mat]] <- processed_data_mat
    print(paste("Processed .mat file for", variable_name_mat))

    #----------------Process .csv files--------------------
    if (grepl("ADMeta", ad_comb)) {
      variable_name_csv <- paste0(dataset, "_", ad_comb, "_csv")
      processed_data_csv <- read_conjfdr_loci_csv(dataset, ad_comb, pleiofdr_results_dir)

      #----------------Extract matching rows from processed_data_mat based on CHR and BP in processed_data_csv--------------------
      extracted_rows <- lapply(1:nrow(processed_data_csv), function(i) {
        row <- processed_data_csv[i, ]
        CHR <- row$chrnum
        BP <- row$chrpos
        idx <- which(data_list[[variable_name_mat]]$CHR == CHR & data_list[[variable_name_mat]]$BP == BP)
        if (length(idx) > 0) {
          return(data_list[[variable_name_mat]][idx, ])
        }
        return(NULL)
      })

      sug_loci_list[[variable_name_csv]] <- do.call(rbind, extracted_rows)

      #----------------Sort sug_loci_list by CHR and BP--------------------
      if (!is.null(sug_loci_list[[variable_name_csv]])) {
        sug_loci_list[[variable_name_csv]] <- sug_loci_list[[variable_name_csv]][order(sug_loci_list[[variable_name_csv]]$CHR, sug_loci_list[[variable_name_csv]]$BP), ]
      }

      print(paste("Processed .csv file for", variable_name_csv))
    }
  }
}

#----------------Save Processed Data--------------------
save(data_list, file = file.path(output_dir, "data_list.RData"))
save(sug_loci_list, file = file.path(output_dir, "sug_loci_list.RData"))

#----------------Get FDR Matrix--------------------
results_df <- get_sug_loci_fdr(sug_loci_list, data_list)

#----------------Save FDR Matrix--------------------
fwrite(results_df, file.path(output_dir, "sug_fdrmat_results.csv"), row.names = FALSE)

#----------------Add Allele Information--------------------
sug_loci_results <- as.data.frame(unique(results_df))
ref <- fread(ref_file)
sug_loci_results <- merge(sug_loci_results, ref[, .(SNP, A1, A2)],
                          by = c("SNP"),
                          all.x = TRUE)
fwrite(sug_loci_results, file.path(output_dir, "sug_fdrmat_results_withallele.csv"), row.names = FALSE)

#----------------LiftOver to Build 38--------------------
sug_loci_results$rs37 <- sug_loci_results$SNP

output_file <- file.path(output_dir, "sug_fdrmat_results")
sug_loci_results <- create_bed_file(sug_loci_results, output_file)
sug_loci_results <- perform_liftOver(sug_loci_results, liftover_dir, output_file)

sug_loci_results$SNP <- sug_loci_results$rs37
sug_loci_results <- sug_loci_results %>% dplyr::select(-A1, -A2, -rs37)
fwrite(sug_loci_results, file.path(output_dir, "sug_fdrmat_results_hg38.csv"), row.names = FALSE)

#----------------Group Loci by LD--------------------
results <- get_ld_groups(sug_loci_results,
                         plink_exe = plink_exe,
                         keep_file = ld_keep_file,
                         bfile_prefix = ld_bfile_prefix,
                         work_dir = ld_work_dir,
                         ld_threshold = ld_group_threshold)
fwrite(results, file.path(output_dir, "sug_fdrmat_results_hg38_grouped.csv"), row.names = FALSE)

#----------------Filter Results with p-value < 0.05--------------------
sug_loci_results <- as.data.table(results)

cols_to_check <- names(sug_loci_results)[which(names(sug_loci_results) == paste0(traits$label[1], "_ADMeta_female_mat")):
                                         which(names(sug_loci_results) == paste0(traits$label[nrow(traits)], "_ADMeta_all_mat"))]

long_data <- melt(sug_loci_results,
                  id.vars = c("SNP", "CHR", "BP", "Group"),
                  measure.vars = cols_to_check,
                  variable.name = "col_name",
                  value.name = "p_value")

filtered_data <- long_data[with(long_data, p_value < 0.05), ]
filtered_data <- filtered_data %>% arrange(CHR, BP)
filtered_data$rs37 <- filtered_data$SNP

fwrite(filtered_data, file.path(output_dir, "sug_fdrmat_results_hg38_grouped_filtered.csv"), row.names = FALSE)
