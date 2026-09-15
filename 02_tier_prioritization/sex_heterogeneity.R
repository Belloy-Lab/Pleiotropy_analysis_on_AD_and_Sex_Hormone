# =============================================================================
# 02_tier_prioritization / sex_heterogeneity.R
#
# Tier-2 component (AD sex heterogeneity). Starting from the LD-grouped conjFDR
# locus table, this script writes:
#
#   sex_het_p.tsv         AD sex-heterogeneity P value (SEX_HET_P) of the index
#                         variant, written into the trait columns of the locus
#                         table; 1 when the variant is not present in the
#                         sex-heterogeneity GWAS
#
#   sex_het_p_1.5fold.tsv same table, with the value 0.04 written for variants
#                         whose female/male AD effect-size ratio exceeds 1.5 in
#                         either direction (the effect-size component of Tier 2)
#
# Column layout expected for the grouped locus table (output of ld_grouping.R):
# the trait columns come first (one column per trait x AD stratum combination),
# followed by the identifier columns CHR, BP, SNP and Group.
# =============================================================================

# =========================
# USER CONFIGURATION
# =========================
ad_female_file    <- "path/to/ad/ad_female_gwas.txt.gz"             # columns: SNP, BETA, ...
ad_male_file      <- "path/to/ad/ad_male_gwas.txt.gz"               # columns: SNP, BETA, ...
sex_het_file      <- "path/to/ad/ad_sex_heterogeneity_gwas.txt.gz"  # columns: SNP, SEX_HET_P
grouped_loci_file <- "path/to/output/sug_fdrmat_results_hg38_grouped.csv"
output_dir        <- "path/to/output/"

# Number of trait columns at the start of the grouped locus table
# (number of traits x number of AD strata; 11 traits x 2 strata = 22 for the manuscript)
n_trait_columns <- 22

# Values used by the original analysis
ratio_threshold    <- 1.5    # female/male AD effect-size ratio flagging the criterion
ratio_flag_value   <- 0.04   # written when the ratio threshold is exceeded
missing_p_value    <- 1      # written when the variant is absent from the input GWAS
# =========================

library(data.table)
library(dplyr)

#------------------Load AD GWAS--------------------
AD_female_meta <- fread(ad_female_file)
AD_male_meta <- fread(ad_male_file)

AD_female_meta <- as.data.frame(AD_female_meta)
AD_male_meta <- as.data.frame(AD_male_meta)

#------------------Load grouped conjFDR locus table--------------------
sug_loci_results <- fread(grouped_loci_file)
sug_loci_results <- as.data.frame(sug_loci_results)
sug_loci_results <- unique(sug_loci_results)

#------------------Get sex het p-values for each SNP--------------------
sex_het_results <- fread(sex_het_file)

#------------------copy to avoid overwriting--------------------
sex_het_p_result <- copy(sug_loci_results)

#------------------Get sex het p-values for each SNP--------------------
for (i in 1:nrow(sex_het_p_result)) {
  snp <- sex_het_p_result$SNP[i]

  snp_sex_het_p <- sex_het_results %>% dplyr::filter(SNP == snp) %>% dplyr::pull(SEX_HET_P)

  if (length(snp_sex_het_p) == 0) {
    sex_het_p_result[i, 1:n_trait_columns] <- missing_p_value
    next
  }

  sex_het_p_result[i, 1:n_trait_columns] <- snp_sex_het_p
}

#------------------Save results--------------------
fwrite(sex_het_p_result, file.path(output_dir, "sex_het_p.tsv"), quote=F, row.names=F, col.names=T, sep='\t')

#------------------Flag a female/male effect-size ratio > 1.5--------------------
for (i in 1:nrow(sex_het_p_result)) {
  snp <- sex_het_p_result$SNP[i]

  snp_sex_beta_female <- AD_female_meta %>% dplyr::filter(SNP == snp) %>% dplyr::pull(BETA)
  snp_sex_beta_male <- AD_male_meta %>% dplyr::filter(SNP == snp) %>% dplyr::pull(BETA)

  if (length(snp_sex_beta_female) == 0 | length(snp_sex_beta_male)==0) {
    sex_het_p_result[i, 1:n_trait_columns] <- missing_p_value
    next
  }

  beta_ratio <- snp_sex_beta_female / snp_sex_beta_male

  if (beta_ratio > ratio_threshold) {
    sex_het_p_result[i, 1:n_trait_columns] <- ratio_flag_value
    next
  }

  beta_ratio <- snp_sex_beta_male / snp_sex_beta_female

  if (beta_ratio > ratio_threshold) {
    sex_het_p_result[i, 1:n_trait_columns] <- ratio_flag_value
    next
  }
}

#------------------Save results--------------------
fwrite(sex_het_p_result, file.path(output_dir, "sex_het_p_1.5fold.tsv"), quote=F, row.names=F, col.names=T, sep='\t')
