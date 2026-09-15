# =============================================================================
# 02_tier_prioritization / coloc_abf.R
#
# AD x sex hormone-related trait colocalization with coloc.abf.
# For every locus/trait pair with conjFDR < 0.05 the script extracts a +/-
# 1 Mb window around the index variant from the AD GWAS and from the hormone
# trait GWAS, harmonises the two datasets and computes the colocalization
# posterior probability (PP4) with coloc.abf (default priors).
#
# Required columns
#   AD GWAS       : CHR, BP, SNP, ALLELE1, ALLELE0, A1FREQ, N_incl, BETA, SE, P
#   hormone trait : CHR, BP, rs37, A1, A2, EAF, BETA, SE, P, N
# =============================================================================

# =========================
# USER CONFIGURATION
# =========================
grouped_loci_file <- "path/to/output/sug_fdrmat_results_hg38_grouped.csv"

ad_female_file <- "path/to/ad/ad_female_gwas.txt.gz"
ad_male_file   <- "path/to/ad/ad_male_gwas.txt.gz"
ad_all_file    <- "path/to/ad/ad_all_gwas.txt.gz"

# One entry per sex hormone-related trait that is present in the grouped locus
# table. The names are the trait labels used in the pleioFDR results (see
# ld_grouping.R) and are matched against the column names of the locus table.
hormone_files <- list(
  Meno   = "path/to/traits/age_at_menopause_hg38.txt",
  Mena   = "path/to/traits/age_at_menarche_hg38.txt",
  voice  = "path/to/traits/age_at_voice_breaking_hg38.txt",
  BTF    = "path/to/traits/bioavailable_testosterone_female_hg38.txt",
  BTM    = "path/to/traits/bioavailable_testosterone_male_hg38.txt",
  TTF    = "path/to/traits/total_testosterone_female_hg38.txt",
  TTM    = "path/to/traits/total_testosterone_male_hg38.txt",
  SBFBMI = "path/to/traits/shbg_female_bmi_adjusted_hg38.txt",
  SBMBMI = "path/to/traits/shbg_male_bmi_adjusted_hg38.txt",
  SBF    = "path/to/traits/shbg_female_hg38.txt",
  SBM    = "path/to/traits/shbg_male_hg38.txt"
)

output_file <- "path/to/output/sug_coloc_PP4_results.csv"
# =========================

library(data.table)
library(dplyr)
library(coloc)

source("functions_coloc_abf.R")

#------------------Load Grouped Loci Results--------------------
sug_loci_results <- fread(grouped_loci_file)

#------------------Load AD Meta-Analysis Data--------------------
AD_Meta_female <- as.data.frame(fread(ad_female_file))
AD_Meta_male   <- as.data.frame(fread(ad_male_file))
AD_Meta_all    <- as.data.frame(fread(ad_all_file))

#------------------Rename Columns--------------------
ad_column_map_from <- c("ALLELE1", "ALLELE0", "A1FREQ", "N_incl", "SNP")
ad_column_map_to   <- c("A1", "A2", "EAF", "N", "rs38")

setnames(AD_Meta_female, ad_column_map_from, ad_column_map_to)
setnames(AD_Meta_male,   ad_column_map_from, ad_column_map_to)
setnames(AD_Meta_all,    ad_column_map_from, ad_column_map_to)

AD_Meta_female$SNP <- AD_Meta_female$rs38
AD_Meta_male$SNP   <- AD_Meta_male$rs38
AD_Meta_all$SNP    <- AD_Meta_all$rs38

#------------------Load Hormone Data--------------------
new_horm_data_list <- list()
for (hormone_label in names(hormone_files)) {
  new_horm_data_list[[hormone_label]] <- fread(hormone_files[[hormone_label]])
  new_horm_data_list[[hormone_label]]$SNP <- new_horm_data_list[[hormone_label]]$rs37
}

#------------------New Data Frame for Storing Results--------------------
sug_loci_results <- as.data.frame(sug_loci_results)
PP4_results <- sug_loci_results %>% copy()

#------------------Run Coloc Analysis--------------------
PP4_results <- run_coloc_analysis(sug_loci_results, AD_Meta_female, AD_Meta_male, AD_Meta_all, new_horm_data_list, PP4_results)
fwrite(PP4_results, output_file, row.names = FALSE)
