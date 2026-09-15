# Publication version of the xQTL colocalization workflow
# supporting "Shared Risk Genes and Causal Relationships across
# Sex Hormone Related Traits and Alzheimer's Disease".
#
# Bulk colocalization driver (coloc.susie) for one QTL dataset identified by
# study / tissue / qtl_type. The script is sourced by xqtl_coloc_susie.R and expects
# the following objects in the calling environment:
#
#   qtl_file         formatted QTL dataset (one row per variant x molecular trait)
#   tier1_loci_file  Tier-1 locus table with columns chrom, bp_38_start, bp_38_end,
#                    bp_38_med, locus_index, stratum, discovery
#   study, tissue, qtl_type, output_dir
#   ref_seq          GENCODE annotation as a data.table
#   AD strata        Female_rb, Male_rb, NS_rb, each holding the case
#                    proportion in a column named `s`
#   plink_bin, ld_reference_prefix, temp_dir   (see xqtl_coloc_susie.R)
#
# The AD strata retained for this manuscript are Female, Male and NonStrat.
# Supported qtl_type values are eQTL, sQTL and pQTL. The mQTL, caQTL and haQTL
# implementations in functions/ require dataset-specific harmonisation and are not
# dispatched by this driver (see README).
#
# Note: this driver uses an allele frequency agreement threshold of 0.5 between AD and
# the QTL dataset, whereas bulk_colocalization_abf.R uses 0.1. The difference is
# original and is preserved here.

library(data.table)
library(dplyr)
library(stringr)

#------------------Read the formatted QTL dataset and the Tier-1 locus table----------
qtl_data <- fread(as.character(qtl_file))
all_suggestive_snps <- fread(as.character(tier1_loci_file))

final_results <- data.frame()
colocalization_results <- data.frame()

for (i in 1:nrow(all_suggestive_snps)) {

  # Locus definition from the Tier-1 input table
  chrom <- as.integer(all_suggestive_snps[i, chrom])
  bp_38_start <- as.double(all_suggestive_snps[i, bp_38_start])
  bp_38_end <- as.double(all_suggestive_snps[i, bp_38_end])
  bp_38_med <- as.double(all_suggestive_snps[i, bp_38_med])
  locus_index <- all_suggestive_snps[i, locus_index]
  stratum <- all_suggestive_snps[i, stratum]
  discovery <- all_suggestive_snps[i, discovery]

  message("Starting data merge of ", locus_index, " in ", study, " ", tissue, " ", qtl_type)

  #------------------Select the AD stratum--------------------
  if (stratum == "Female") {
    ad_df <- Female_rb
  } else if (stratum == "Male") {
    ad_df <- Male_rb
  } else if (stratum == "NonStrat") {
    ad_df <- NS_rb
  } else {
    stop("Unknown AD stratum '", stratum, "': expected Female, Male or NonStrat")
  }

  # AD GWAS restricted to the chromosome of the locus
  ad_df2 <- ad_df %>%
    filter(CHR == chrom)

  # GENCODE genes on this chromosome, one record per gene (longest transcript)
  ref_seq_2 <- ref_seq %>%
    filter(type == "gene", seqnames == paste0("chr", chrom)) %>%
    group_by(gene_name) %>%
    filter(width == max(width)) %>%
    ungroup() %>%
    arrange(start) %>%
    dplyr::select(start, end, width, gene_name, gene_id, transcript_id, gene_type) %>%
    mutate(gene_id = str_extract(gene_id, "^[^.]+"))
  setDT(ref_seq_2)

  #------------------QTL dataset restricted to the locus chromosome--------------------
  df_2 <- qtl_data %>%
    filter(CHR == chrom) %>%
    arrange(BP)

  df_3 <- merge(df_2, ref_seq_2, all.x = TRUE, by = "gene_id") %>%
    arrange(BP)

  if (study == "MetaBrain" | study == "eQTLgen") {
    df_3 <- df_3 %>%
      mutate(gene_name = GeneSymbol)
  }

  #------------------Harmonise AD GWAS and QTL dataset on the variant identifier--------------------
  if (study == "kosoy" | study == "eQTLgen" | study == "SingleBrain") {

    merged_12 <- merge(df_3, ad_df2, by.x = "id12", by.y = "posID")
    merged_21 <- merge(df_3, ad_df2, by.x = "id21", by.y = "posID") %>%
      mutate(beta = beta * -1)

    merged_updated <- rbind(merged_12, merged_21) %>%
      arrange(CHR.x, BP.x)

    merged_updated <- merged_updated %>%
      dplyr::rename(CHR = CHR.x, BP = BP.x, AD_P = P, AD_EAF = A1FREQ, AD_BETA = BETA,
                    AD_SE = SE, AD_N = N_incl, QTL_P = pvalue, QTL_BETA = beta,
                    QTL_SE = se, QTL_N = N) %>%
      mutate(QTL_EAF = AD_EAF) %>%
      dplyr::select(CHR, BP, SNP, ALLELE1, ALLELE0, AD_EAF, AD_BETA, AD_SE, Z, AD_P, AD_N,
                    QTL_EAF, QTL_BETA, QTL_SE, QTL_P, QTL_N, gene_name,
                    molecular_trait_object_id, s)

  } else {

    merged_12 <- merge(df_3, ad_df2, by.x = "id12", by.y = "posID")
    merged_21 <- merge(df_3, ad_df2, by.x = "id21", by.y = "posID") %>%
      mutate(beta = beta * -1,
             EAF = 1 - EAF)

    merged_updated <- rbind(merged_12, merged_21) %>%
      arrange(CHR.x, BP.x)

    merged_updated <- merged_updated %>%
      filter(abs(A1FREQ - EAF) <= 0.5) %>%  # effect allele frequency agreement between AD and QTL
      dplyr::rename(CHR = CHR.x, BP = BP.x, AD_P = P, AD_EAF = A1FREQ, AD_BETA = BETA,
                    AD_SE = SE, AD_N = N_incl, QTL_P = pvalue, QTL_EAF = EAF,
                    QTL_BETA = beta, QTL_SE = se, QTL_N = N) %>%
      dplyr::select(CHR, BP, SNP, ALLELE1, ALLELE0, AD_EAF, AD_BETA, AD_SE, Z, AD_P, AD_N,
                    QTL_EAF, QTL_BETA, QTL_SE, QTL_P, QTL_N, gene_name,
                    molecular_trait_object_id, s)
  }

  merged <- merged_updated

  name_dset <- paste(study, tissue, "chr", chrom, qtl_type, "RB", stratum, sep = "_")

  rm(merged_12, merged_21, ad_df, ad_df2, df_3, ref_seq_2)
  gc()

  #------------------Dispatch to the QTL-type specific colocalization function--------------------
  if (qtl_type == "eQTL") {
    source("functions/eqtl_susie.R")
    result <- eQTL_colocalization(merged = merged_updated,
                                  all_suggestive_snps = all_suggestive_snps,
                                  name_dset = name_dset)
  } else if (qtl_type == "sQTL") {
    source("functions/sqtl_susie.R")
    result <- sQTL_colocalization(merged = merged_updated,
                                  all_suggestive_snps = all_suggestive_snps,
                                  name_dset = name_dset)
  } else if (qtl_type == "pQTL") {
    source("functions/pqtl_susie.R")
    result <- pQTL_colocalization(merged = merged_updated,
                                  all_suggestive_snps = all_suggestive_snps,
                                  name_dset = name_dset)
  } else {
    stop("qtl_type '", qtl_type, "' is not dispatched by this driver. The mQTL, caQTL ",
         "and haQTL implementations in functions/ require dataset-specific harmonisation ",
         "(see README).")
  }

  if (is.null(result)) {
    message("No results for: ", locus_index, " in tissue: ", tissue, " ", qtl_type)
    next
  }

  #------------------Accumulate results--------------------
  colocalization_results <- rbind(colocalization_results, result)

  if (nrow(colocalization_results) > 0) {
    final_results <- rbind(final_results, colocalization_results)
  }
}

#------------------Arrange the final results and write one results table--------------------
final_results <- final_results %>%
  arrange(desc(PP4)) %>%
  group_by(locus, gene_name, PP3, PP4) %>%
  distinct(PP4, .keep_all = TRUE)
final_results <- as.data.frame(final_results)

fwrite(final_results, file.path(output_dir, paste(study, tissue, qtl_type, "susie_hg38.csv", sep = "_")))

message("Completed QTL colocalization SuSiE analysis for: ", study, " ", tissue)
