# Publication version of the xQTL colocalization workflow
# supporting "Shared Risk Genes and Causal Relationships across
# Sex Hormone Related Traits and Alzheimer's Disease".
#
# haQTL (histone acetylation) colocalization with coloc.susie.
#
# Inputs
#   merged              harmonised AD x haQTL table for one locus with columns
#                       CHR, BP, SNP, ALLE1, ALLE0, AD_EAF, AD_BETA, AD_SE, AD_P,
#                       AD_N, QTL_EAF, QTL_BETA, QTL_SE, QTL_P, QTL_N, peak, peakPos, s
#   all_suggestive_snps Tier-1 locus table (accepted for interface compatibility; unused)
#   name_dset           dataset label written to the `dset` column of the result table
#
# Per-locus variables are taken from the calling environment, as in the original
# pipeline: chrom, bp_38_start, bp_38_end, bp_38_med, locus_index, study, tissue,
# discovery.
#
# LD reference configuration, defined by the entry-point script:
#   plink_bin            PLINK executable
#   ld_reference_prefix  PLINK binary prefix of the LD reference panel; the chromosome
#                        number is appended (e.g. "path/to/LD_reference/TOPMed_chr" + 1)
#   temp_dir             directory used for the per-locus SNP list and PLINK LD files
#
# Returns one row per histone-acetylation peak with PP0-PP4 from the credible-set pair
# with the largest PP4.

haQTL_colocalization <- function(merged, all_suggestive_snps, name_dset) {

  bp <- as.integer(bp_38_med)
  bp_s <- as.integer(bp_38_start)
  bp_e <- as.integer(bp_38_end)
  chrom <- as.integer(chrom)
  range <- 1e6

  if (discovery != "pleiotropy") {
    colocalizing_genes <- merged %>%
      filter(BP >= bp - range & BP <= bp + range) %>%
      distinct(peak, .keep_all = TRUE) %>%
      arrange(peak)
  } else {
    colocalizing_genes <- merged %>%
      filter(BP >= bp_s - range & BP <= bp_e + range) %>%
      distinct(peak, .keep_all = TRUE) %>%
      arrange(peak)
  }

  if (nrow(colocalizing_genes) == 0) {
    message("No peaks in the analysis window for ", name_dset)
    return(NULL)
  }

  colocalization_results <- data.frame()

  for (gene_row in 1:nrow(colocalizing_genes)) {

    gene_n <- colocalizing_genes$peak[gene_row]
    gene_pos <- colocalizing_genes$peakPos[gene_row]

    if (discovery != "pleiotropy") {
      merged_c <- merged %>%
        mutate(varbeta_ad = AD_SE^2,
               varbeta_t2 = QTL_SE^2) %>%
        filter(BP >= (bp - range) & BP <= (bp + range)) %>%
        filter(varbeta_ad != 0, varbeta_t2 != 0, !is.na(SNP)) %>%
        filter(peak == gene_n,
               QTL_P > 0) %>%
        distinct(SNP, .keep_all = TRUE)
    } else {
      merged_c <- merged %>%
        mutate(varbeta_ad = AD_SE^2,
               varbeta_t2 = QTL_SE^2) %>%
        filter(BP >= (bp_s - range) & BP <= (bp_e + range)) %>%
        filter(varbeta_ad != 0, varbeta_t2 != 0, !is.na(SNP)) %>%
        filter(peak == gene_n,
               QTL_P > 0) %>%
        distinct(SNP, .keep_all = TRUE)
    }

    if (nrow(merged_c) == 0) {
      message("No SNPs for ", gene_n, " in the merged file")
      next
    }

    snps_for_ld <- list(merged_c$SNP)
    snp_file <- file.path(temp_dir, paste0(study, "_", tissue, "_haQTL_", discovery, "_", locus_index, ".txt"))
    fwrite(snps_for_ld, snp_file, sep = "\t", col.names = FALSE)

    file_ld <- paste("LD_matrix", study, tissue, "haQTL", discovery, locus_index, sep = "_")
    ld_prefix <- file.path(temp_dir, file_ld)
    ld_data_prefix <- paste0(ld_prefix, ".data")

    command <- paste(plink_bin,
                     "--bfile", paste0(ld_reference_prefix, chrom),
                     "--allow-no-sex",
                     "--extract", snp_file,
                     "--keep-allele-order",
                     "--make-bed --out", ld_data_prefix,
                     sep = " ")
    system(command)

    command <- paste(plink_bin,
                     "--bfile", ld_data_prefix,
                     "--allow-no-sex",
                     "--keep-allele-order",
                     "--r --matrix",
                     "--out", ld_prefix,
                     sep = " ")
    system(command)

    bim <- fread(paste0(ld_data_prefix, ".bim"))
    snps <- unlist(bim$V2)

    LD <- fread(paste0(ld_prefix, ".ld"))
    LD <- as.matrix(LD)
    colnames(LD) <- snps
    rownames(LD) <- snps

    file.remove(paste0(ld_data_prefix, ".bed"))
    file.remove(paste0(ld_data_prefix, ".bim"))
    file.remove(paste0(ld_data_prefix, ".fam"))
    file.remove(paste0(ld_prefix, ".ld"))
    file.remove(paste0(ld_prefix, ".log"))
    file.remove(paste0(ld_prefix, ".nosex"))
    file.remove(paste0(ld_data_prefix, ".log"))
    file.remove(paste0(ld_data_prefix, ".nosex"))
    file.remove(snp_file)

    ADl <- as.list(merged_c[, c("AD_BETA", "varbeta_ad", "SNP", "BP", "AD_EAF")])
    names(ADl)[1:5] <- c("beta", "varbeta", "snp", "position", "MAF")
    names(ADl$beta) <- ADl$snp
    names(ADl$varbeta) <- ADl$snp
    ADl$N <- median(merged$AD_N, na.rm = TRUE)

    matching_snps_ADl <- ADl$snp %in% snps
    ADl$beta <- ADl$beta[matching_snps_ADl]
    ADl$varbeta <- ADl$varbeta[matching_snps_ADl]
    ADl$snp <- ADl$snp[matching_snps_ADl]
    ADl$position <- ADl$position[matching_snps_ADl]
    ADl$MAF <- ADl$MAF[matching_snps_ADl]
    ADl$LD <- LD
    ADl$type <- "cc"
    cc <- ADl

    common_snps <- intersect(cc$snp, colnames(cc$LD))
    cc_index <- match(common_snps, cc$snp)
    ld_index <- match(common_snps, colnames(cc$LD))

    cc$beta <- cc$beta[cc_index]
    cc$varbeta <- cc$varbeta[cc_index]
    cc$snp <- cc$snp[cc_index]
    cc$position <- cc$position[cc_index]
    cc$MAF <- cc$MAF[cc_index]
    cc$LD <- cc$LD[ld_index, ld_index]

    traitl <- as.list(merged_c[, c("QTL_BETA", "varbeta_t2", "SNP", "BP", "QTL_EAF")])
    names(traitl)[1:5] <- c("beta", "varbeta", "snp", "position", "MAF")
    names(traitl$beta) <- traitl$snp
    names(traitl$varbeta) <- traitl$snp
    traitl$N <- round(median(merged$QTL_N, na.rm = TRUE))

    matching_snps_traitl <- traitl$snp %in% snps
    traitl$beta <- traitl$beta[matching_snps_traitl]
    traitl$varbeta <- traitl$varbeta[matching_snps_traitl]
    traitl$snp <- traitl$snp[matching_snps_traitl]
    traitl$position <- traitl$position[matching_snps_traitl]
    traitl$MAF <- traitl$MAF[matching_snps_traitl]
    traitl$LD <- LD
    traitl$type <- "quant"
    eq <- traitl

    common_snps <- intersect(eq$snp, colnames(eq$LD))
    eq_index <- match(common_snps, eq$snp)
    ld_index <- match(common_snps, colnames(eq$LD))

    eq$beta <- eq$beta[eq_index]
    eq$varbeta <- eq$varbeta[eq_index]
    eq$snp <- eq$snp[eq_index]
    eq$position <- eq$position[eq_index]
    eq$MAF <- eq$MAF[eq_index]
    eq$LD <- eq$LD[ld_index, ld_index]

    if (length(eq$LD) == 1 | length(cc$LD) == 1) {
      next
    }

    coverage_values <- seq(0.95, 0.05, by = -0.1)
    cs_length_S1 <- 0
    cs_length_S2 <- 0
    for (coverage in coverage_values) {
      tryCatch({
        S1 <- runsusie(cc, coverage = coverage, max_iter = 10000)
        cs_length_S1 <- length(S1$sets$cs)

        S2 <- runsusie(eq, coverage = coverage, max_iter = 10000)
        cs_length_S2 <- length(S2$sets$cs)

        if (cs_length_S1 > 0 & cs_length_S2 > 0) {
          break
        }
      }, error = function(e) {
        message("Error encountered: ", e)
      })
    }

    gc()

    if (!exists("S1") || !exists("S2")) {
      message("S1 or S2 does not exist for: ", gene_n)
      next
    }

    if (class(S1) != "susie" || class(S2) != "susie") {
      message("S1 or S2 failed to converge for: ", gene_n)
      next
    }

    out <- coloc.susie(S1, S2)

    coloc_PP0 <- 0
    coloc_PP1 <- 0
    coloc_PP2 <- 0
    coloc_PP3 <- 0
    coloc_PP4 <- 0
    coloc_hit1 <- 0
    coloc_hit2 <- 0

    if (is.null(out$summary$PP.H4.abf)) {
      message("No coloc results for: ", gene_n)
    } else {
      index_row <- which.max(out$summary$PP.H4.abf)
      coloc_PP0 <- out$summary$PP.H0.abf[index_row]
      coloc_PP1 <- out$summary$PP.H1.abf[index_row]
      coloc_PP2 <- out$summary$PP.H2.abf[index_row]
      coloc_PP3 <- out$summary$PP.H3.abf[index_row]
      coloc_PP4 <- out$summary$PP.H4.abf[index_row]

      if (is.na(coloc_PP4)) {
        coloc_PP4 <- 0
      }

      coloc_hit1 <- out$summary$hit1[index_row]
      coloc_hit2 <- out$summary$hit2[index_row]
    }

    colocalization_row <- data.frame(
      CHROM = chrom,
      BP = bp,
      dset = name_dset,
      locus = locus_index,
      gene_name = gene_n,
      molecular_trait_id = "",
      PP0 = coloc_PP0,
      PP1 = coloc_PP1,
      PP2 = coloc_PP2,
      PP3 = coloc_PP3,
      PP4 = coloc_PP4,
      qtl_type = "haQTL",
      method = "susie",
      hit1 = coloc_hit1,
      hit2 = coloc_hit2,
      n_snps = nrow(LD),
      discovery = discovery,
      stringsAsFactors = FALSE
    )

    colocalization_results <- rbind(colocalization_results, colocalization_row)
  }

  return(colocalization_results)
}
