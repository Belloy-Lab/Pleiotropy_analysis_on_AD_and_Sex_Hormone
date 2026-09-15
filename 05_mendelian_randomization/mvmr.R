# =============================================================================
# 05_mendelian_randomization / mvmr.R
#
# Multivariable MR (MVMR) of sex hormone-binding globulin (SHBG) and BMI on
# sex-stratified Alzheimer's disease, using the MVMR implementations provided
# by TwoSampleMR.
# =============================================================================

# =========================
# USER CONFIGURATION
# =========================
outcome_rdata   <- "path/to/ad_outcome_data.RData"        # contains `trait_all` (see README)
bmi_input_file  <- "path/to/mvmr_input/bmi_male_mvmr_input.txt"
shbg_input_file <- "path/to/mvmr_input/shbg_male_mvmr_input.txt"
ld_bfile        <- "path/to/plink_reference/EUR"           # PLINK binary prefix (EUR)
output_dir      <- "path/to/output/mvmr/"

mvmr_pval_threshold <- 5e-08
clump_r2            <- 0.001
clump_kb            <- 10000
ad_stratum          <- "AD male"                           # outcome label used for MVMR
# =========================

library(TwoSampleMR)
library(vroom)
library(dplyr)
library(ieugwasr)
library(plinkbinr)
library(data.table)
plink_pathway <- get_plink_exe()

# Load AD dataset
load(outcome_rdata)       # provides `trait_all`
outcome_data <- trait_all
# rename column
colnames(outcome_data) <- c("chr.outcome",
                            "pos.outcome",
                            "SNP",
                            "effect_allele.outcome",
                            "other_allele.outcome",
                            "eaf.outcome",
                            "beta.outcome",
                            "se.outcome",
                            "pval.outcome",
                            "samplesize.outcome",
                            "outcome")
outcome_data$id.outcome <- outcome_data$outcome

# create MVMR IV file
exposure_data_mvmr = mv_extract_exposures_local(
                     c(bmi_input_file,
                       shbg_input_file),
                     sep = "\t",
                     phenotype_col = "exposure",
                     snp_col = "SNP",
                     beta_col = "beta.exposure",
                     se_col = "se.exposure",
                     eaf_col = "eaf.exposure",
                     effect_allele_col = "effect_allele.exposure",
                     other_allele_col = "other_allele.exposure",
                     pval_col = "pval.exposure",
                     id_col = "id.exposure",
                     min_pval = 1e-200,
                     log_pval = FALSE,
                     pval_threshold = mvmr_pval_threshold,
                     plink_bin = plink_pathway,
                     bfile = ld_bfile,
                     clump_r2 = clump_r2,
                     clump_kb = clump_kb,
                     pop = "EUR",
                     harmonise_strictness = 2)

# Filter the AD GWAS dataset for the MVMR analysis
outcome_data_mvmr <- filter(outcome_data, outcome_data$id.outcome == ad_stratum)

# Harmonization for MVMR
mvdat <- mv_harmonise_data(exposure_dat = exposure_data_mvmr,
                           outcome_dat = outcome_data_mvmr,
                           harmonise_strictness = 3)
save(mvdat, file = file.path(output_dir, "mvdat.RData"))

# MVMR
methods <- c("mv_basic", "mv_ivw", "mv_multiple", "mv_residual")

# Run each method and collect results
mvmr_res_list <- lapply(methods, function(method) {
  # Dynamically call the method
  mvmr_res <- do.call(method, list(mvdat, pval_threshold = mvmr_pval_threshold))

  # Extract the result and generate odds ratios
  mvmr_res_df <- mvmr_res[["result"]]
  mvmr_res_df <- generate_odds_ratios(mvmr_res_df)

  # Add additional information columns
  mvmr_res_df$mvmr_function <- method

  return(mvmr_res_df)
})

mvmr_res <- do.call(rbind, mvmr_res_list)

write.csv(mvmr_res, file.path(output_dir, "mvmr_res.csv"), row.names = F)
