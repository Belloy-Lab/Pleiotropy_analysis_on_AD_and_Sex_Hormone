# Publication version of the xQTL colocalization workflow
# supporting "Shared Risk Genes and Causal Relationships across
# Sex Hormone Related Traits and Alzheimer's Disease".
#
# pQTL colocalization with coloc.abf.
#
# Inputs
#   merged              harmonised AD x pQTL table for one locus with columns
#                       CHR, BP, SNP, ALLELE1, ALLELE0, AD_EAF, AD_BETA, AD_SE, AD_P,
#                       AD_N, QTL_EAF, QTL_BETA, QTL_SE, QTL_P, QTL_N, gene_name, s
#   all_suggestive_snps Tier-1 locus table (accepted for interface compatibility; unused)
#   name_dset           dataset label written to the `dset` column of the result table
#
# Per-locus variables are taken from the calling environment, as in the original
# pipeline: chrom, bp_38_start, bp_38_end, bp_38_med, locus_index, study, tissue,
# discovery.
#
# Returns one row per protein (gene) with the posterior probabilities PP0-PP4.

pQTL_colocalization <- function(merged, all_suggestive_snps, name_dset) {

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

    if (discovery != "pleiotropy") {
      merged_c <- merged %>%
        filter(BP >= (bp - range) & BP <= (bp + range)) %>%
        filter(AD_P > 0,
               QTL_P > 0,
               gene_name == gene_n) %>%
        mutate(varbeta_ad = AD_SE^2,
               varbeta_t2 = QTL_SE^2) %>%
        distinct(SNP, .keep_all = TRUE)
    } else {
      merged_c <- merged %>%
        filter(BP >= (bp_s - range) & BP <= (bp_e + range)) %>%
        filter(AD_P > 0,
               QTL_P > 0,
               gene_name == gene_n) %>%
        mutate(varbeta_ad = AD_SE^2,
               varbeta_t2 = QTL_SE^2) %>%
        distinct(SNP, .keep_all = TRUE)
    }

    if (nrow(merged_c) == 0) {
      next
    }

    D1 <- list(
      pvalues = merged_c$AD_P,
      type = "cc",
      snp = merged_c$SNP,
      position = merged_c$BP,
      N = merged_c$AD_N,
      MAF = merged_c$AD_EAF,
      s = merged_c$s[1]
    )
    check_dataset(D1)

    D2 <- list(
      pvalues = merged_c$QTL_P,
      type = "quant",
      snp = merged_c$SNP,
      position = merged_c$BP,
      N = merged_c$QTL_N,
      MAF = merged_c$QTL_EAF
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
      molecular_trait_id = "",
      PP0 = pp_values["PP.H0.abf"],
      PP1 = pp_values["PP.H1.abf"],
      PP2 = pp_values["PP.H2.abf"],
      PP3 = pp_values["PP.H3.abf"],
      PP4 = pp_values["PP.H4.abf"],
      qtl_type = "pQTL",
      method = "abf",
      hit1 = 0,
      hit2 = 0,
      n_snps = n_distinct(merged_c$SNP),
      discovery = discovery,
      stringsAsFactors = FALSE
    )

    colocalization_results <- rbind(colocalization_results, colocalization_row)
  }

  return(colocalization_results)
}
