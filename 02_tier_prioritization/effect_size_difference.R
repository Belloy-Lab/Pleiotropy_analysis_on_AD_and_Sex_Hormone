# =============================================================================
# 02_tier_prioritization / effect_size_difference.R
#
# Tier-2 component (AD effect-size difference between sexes). For a table of
# candidate loci with female and male AD effect sizes (and the corresponding
# effect sizes from the no-UKB replication), this script
#
#   1. computes the female, male and sex-difference effect-size deltas between
#      the two AD GWAS versions
#   2. joins the maximum colocalization PP4 (PP4_max) per variant and keeps the
#      variant with the highest PP4_max per LD group
#   3. joins the AD sex-heterogeneity P value per variant and keeps one variant
#      per LD group
#   4. reports the squared correlation (R2) between the two AD GWAS versions for
#      the female, male and sex-difference effect sizes, overall and for the two
#      per-group selections
#
# The original analysis script annotated these R2 values on scatter plots; plots
# are not part of this repository, so the same R2 values are written to a table.
# Required input columns: SNP, Group, BETA_Female, BETA_Male, noUKB_female_BETA,
# noUKB_male_BETA.
# =============================================================================

# =========================
# USER CONFIGURATION
# =========================
candidate_loci_file <- "path/to/input/tier_check_with_betas.csv"       # see required columns above
coloc_results_file  <- "path/to/output/sug_coloc_PP4_results.csv"      # long format: SNP, pairs, PP4, PP4_max
sex_het_file        <- "path/to/ad/ad_sex_heterogeneity_gwas.txt.gz"   # columns: SNP, SEX_HET_P
output_dir          <- "path/to/output/"
# =========================

library(data.table)
library(dplyr)

#------------------Load Data--------------------
sug_check <- fread(candidate_loci_file)

#------------------compute deltas--------------------
sug_check <- sug_check %>%
  mutate(
    beta_female_delta = BETA_Female - noUKB_female_BETA,
    beta_male_delta   = BETA_Male - noUKB_male_BETA,
    beta_sex_delta    = (BETA_Female - BETA_Male) - (noUKB_female_BETA - noUKB_male_BETA)
  )

#------------------Add PP4 information--------------------
PP4_long_results <- fread(coloc_results_file)

sug_check <- sug_check %>%
  left_join(
    PP4_long_results %>%
      group_by(SNP) %>%
      summarise(PP4_max = max(PP4_max, na.rm = TRUE), .groups = "drop"),
    by = "SNP"
  )

#------------------ Compute best PP4 per group --------------------
sug_check_best_PP4_bygroup <- sug_check %>%
  group_by(Group) %>%
  slice_max(PP4_max, n = 1, with_ties = FALSE) %>%
  ungroup()

#------------------ Add sex het P-values --------------------
sex_het_results <- fread(sex_het_file)

sug_check <- sug_check %>%
  left_join(
    sex_het_results %>%
      group_by(SNP) %>%
      summarise(sex_het_p_max = min(SEX_HET_P, na.rm = TRUE), .groups = "drop"),
    by = "SNP"
  )

#------------------ Compute best sex het per group --------------------
sug_check_best_sex_het_bygroup <- sug_check %>%
  group_by(Group) %>%
  slice_max(sex_het_p_max, n = 1, with_ties = FALSE) %>%
  ungroup()

#------------------ Function to compute R squared --------------------
get_r2 <- function(x, y) {
  cor(x, y, use = "complete.obs")^2
}

#------------------ R2 between the two AD GWAS versions --------------------
r2_summary <- rbind(
  data.frame(subset = "all_candidates", comparison = "BETA_Female",
             r2 = get_r2(sug_check$BETA_Female, sug_check$noUKB_female_BETA)),
  data.frame(subset = "all_candidates", comparison = "BETA_Male",
             r2 = get_r2(sug_check$BETA_Male, sug_check$noUKB_male_BETA)),
  data.frame(subset = "all_candidates", comparison = "BETA_sex_difference",
             r2 = get_r2(sug_check$BETA_Female - sug_check$BETA_Male,
                         sug_check$noUKB_female_BETA - sug_check$noUKB_male_BETA)),
  data.frame(subset = "best_PP4_by_group", comparison = "BETA_Female",
             r2 = get_r2(sug_check_best_PP4_bygroup$BETA_Female, sug_check_best_PP4_bygroup$noUKB_female_BETA)),
  data.frame(subset = "best_PP4_by_group", comparison = "BETA_Male",
             r2 = get_r2(sug_check_best_PP4_bygroup$BETA_Male, sug_check_best_PP4_bygroup$noUKB_male_BETA)),
  data.frame(subset = "best_PP4_by_group", comparison = "BETA_sex_difference",
             r2 = get_r2(sug_check_best_PP4_bygroup$BETA_Female - sug_check_best_PP4_bygroup$BETA_Male,
                         sug_check_best_PP4_bygroup$noUKB_female_BETA - sug_check_best_PP4_bygroup$noUKB_male_BETA)),
  data.frame(subset = "best_sex_het_by_group", comparison = "BETA_Female",
             r2 = get_r2(sug_check_best_sex_het_bygroup$BETA_Female, sug_check_best_sex_het_bygroup$noUKB_female_BETA)),
  data.frame(subset = "best_sex_het_by_group", comparison = "BETA_Male",
             r2 = get_r2(sug_check_best_sex_het_bygroup$BETA_Male, sug_check_best_sex_het_bygroup$noUKB_male_BETA)),
  data.frame(subset = "best_sex_het_by_group", comparison = "BETA_sex_difference",
             r2 = get_r2(sug_check_best_sex_het_bygroup$BETA_Female - sug_check_best_sex_het_bygroup$BETA_Male,
                         sug_check_best_sex_het_bygroup$noUKB_female_BETA - sug_check_best_sex_het_bygroup$noUKB_male_BETA))
)

#---------------------Save results--------------------
fwrite(sug_check, file.path(output_dir, "tier_check_with_delta.csv"))
fwrite(sug_check_best_PP4_bygroup, file.path(output_dir, "tier_check_with_delta_best_PP4_by_group.csv"))
fwrite(sug_check_best_sex_het_bygroup, file.path(output_dir, "tier_check_with_delta_best_sex_het_by_group.csv"))
fwrite(r2_summary, file.path(output_dir, "beta_comparison_r2.csv"))
