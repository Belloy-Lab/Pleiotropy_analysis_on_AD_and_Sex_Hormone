# Publication version of the xQTL colocalization workflow
# supporting "Shared Risk Genes and Causal Relationships across
# Sex Hormone Related Traits and Alzheimer's Disease".
#
# sQTL colocalization with coloc.abf.
#
# Inputs
#   merged              harmonised AD x sQTL table for one locus with columns
#                       CHR, BP, SNP, ALLELE1, ALLELE0, AD_EAF, AD_BETA, AD_SE, AD_P,
#                       AD_N, QTL_EAF, QTL_BETA, QTL_SE, QTL_P, QTL_N, gene_name,
#                       molecular_trait_object_id, s
#   all_suggestive_snps Tier-1 locus table (accepted for interface compatibility; unused)
#   name_dset           dataset label written to the `dset` column of the result table
#
# Per-locus variables are taken from the calling environment, as in the original
# pipeline: chrom, bp_38_start, bp_38_end, bp_38_med, locus_index, study, tissue,
# discovery.
#
# Returns one row per intron cluster (gene) with the posterior probabilities PP0-PP4.

sQTL_colocalization <- function(merged, all_suggestive_snps, name_dset) {

  bp <- as.integer(bp_38_med)
  bp_s <- as.integer(bp_38_start)
  bp_e <- as.integer(bp_38_end)
  chrom <- as.integer(chrom)
  range <- 1e6

  if (discovery != "pleiotropy") {
    colocalizing_genes <- merged %>%
      filter(BP >= bp - range & BP <= bp + range) %>%
      distinct(gene_name, .keep_all = TRUE) %>%
      arrange(gene_name) %>%
      filter(!is.na(gene_name))
  } else {
    colocalizing_genes <- merged %>%
      filter(BP >= bp_s - range & BP <= bp_e + range) %>%
      distinct(gene_name, .keep_all = TRUE) %>%
      arrange(gene_name) %>%
      filter(!is.na(gene_name))
  }

  if (nrow(colocalizing_genes) == 0) {
    message("No genes in the analysis window for ", name_dset)
    return(NULL)
  }

  colocalization_results <- data.frame()

  for (gene_row in 1:nrow(colocalizing_genes)) {

    gene_n <- colocalizing_genes$gene_name[gene_row]
    molecular_id <- colocalizing_genes$molecular_trait_object_id[gene_row]

    if (discovery != "pleiotropy") {
      ad_sub <- merged %>%
        mutate(varbeta_ad = AD_SE^2) %>%
        filter(BP >= (bp - range) & BP <= (bp + range)) %>%
        filter(varbeta_ad != 0, !is.na(SNP)) %>%
        filter(gene_name == gene_n, AD_P > 0) %>%
        distinct(SNP, .keep_all = TRUE)

      qtl_sub <- merged %>%
        mutate(varbeta_qtl = QTL_SE^2) %>%
        filter(BP >= (bp - range) & BP <= (bp + range)) %>%
        filter(varbeta_qtl != 0, !is.na(SNP)) %>%
        filter(gene_name == gene_n, QTL_P > 0) %>%
        distinct(SNP, .keep_all = TRUE)
    } else {
      ad_sub <- merged %>%
        mutate(varbeta_ad = AD_SE^2) %>%
        filter(BP >= (bp_s - range) & BP <= (bp_e + range)) %>%
        filter(varbeta_ad != 0, !is.na(SNP)) %>%
        filter(gene_name == gene_n, AD_P > 0) %>%
        distinct(SNP, .keep_all = TRUE)

      qtl_sub <- merged %>%
        mutate(varbeta_qtl = QTL_SE^2) %>%
        filter(BP >= (bp_s - range) & BP <= (bp_e + range)) %>%
        filter(varbeta_qtl != 0, !is.na(SNP)) %>%
        filter(gene_name == gene_n, QTL_P > 0) %>%
        distinct(SNP, .keep_all = TRUE)
    }

    if (nrow(ad_sub) == 0 || nrow(qtl_sub) == 0) {
      next
    }

    D1 <- list(
      pvalues = ad_sub$AD_P,
      type = "cc",
      snp = ad_sub$SNP,
      position = ad_sub$BP,
      N = ad_sub$AD_N,
      MAF = ad_sub$AD_EAF,
      s = ad_sub$s[1]
    )
    check_dataset(D1)

    D2 <- list(
      pvalues = qtl_sub$QTL_P,
      type = "quant",
      snp = qtl_sub$SNP,
      position = qtl_sub$BP,
      N = qtl_sub$QTL_N,
      MAF = qtl_sub$QTL_EAF
    )
    check_dataset(D2)

    coloc_results <- coloc.abf(dataset1 = D1, dataset2 = D2)

    pp_values <- coloc_results$summary[c("PP.H0.abf", "PP.H1.abf", "PP.H2.abf",
                                         "PP.H3.abf", "PP.H4.abf")]
    pp_values[is.nan(pp_values)] <- 0

    colocalization_row <- data.frame(
      CHROM = chrom,
      BP = bp,
      dset = name_dset,
      locus = locus_index,
      gene_name = gene_n,
      molecular_trait_id = molecular_id,
      PP0 = pp_values["PP.H0.abf"],
      PP1 = pp_values["PP.H1.abf"],
      PP2 = pp_values["PP.H2.abf"],
      PP3 = pp_values["PP.H3.abf"],
      PP4 = pp_values["PP.H4.abf"],
      qtl_type = "sQTL",
      method = "abf",
      hit1 = 0,
      hit2 = 0,
      n_snps = n_distinct(ad_sub$SNP),
      discovery = discovery,
      stringsAsFactors = FALSE
    )

    colocalization_results <- rbind(colocalization_results, colocalization_row)
  }

  return(colocalization_results)
}
