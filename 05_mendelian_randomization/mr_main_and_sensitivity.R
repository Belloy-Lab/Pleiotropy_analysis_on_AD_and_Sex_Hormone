# =============================================================================
# 05_mendelian_randomization / mr_main_and_sensitivity.R
#
# Two-sample MR of sex hormone-related traits on Alzheimer's disease with
# sensitivity analyses:
#   * instrument strength (F-statistic)
#   * RadialMR outlier removal before the main analysis
#   * IVW as the primary estimator
#   * Cochran's Q and MR-Egger intercept
#   * MR-Egger, weighted median and weighted mode for FDR-significant results
#   * leave-one-out estimates
#   * MR-PRESSO (global test, outlier test, distortion test)
# =============================================================================

# =========================
# USER CONFIGURATION
# =========================
# Instruments after clumping: RData containing `exposure_data`
instruments_rdata <- "path/to/instruments_after_clumping.RData"
# AD outcome data: RData containing `trait_all` (see README for column order)
outcome_rdata     <- "path/to/ad_outcome_data.RData"
output_dir        <- "path/to/output/mr/"

fdr_method        <- "fdr"     # p.adjust method used for the FDR correction
fdr_threshold     <- 0.05
mrpresso_nb_dist  <- 1000
mrpresso_signif   <- 0.05
# =========================

library(TwoSampleMR)
library(dplyr)
library(MRPRESSO)
library(tidyr)
library(RadialMR)
library(stringr)

# load instruments and outcome data
load(instruments_rdata)   # provides `exposure_data`
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

# add F-statistics
exposure_data <- add_rsq(exposure_data)
exposure_data$F_statistics <- (exposure_data$samplesize.exposure-1-1)*(exposure_data$rsq.exposure/(1-exposure_data$rsq.exposure))

###########################################################
###################### MR analysis ########################
###########################################################

# Harmonize data
H_data <- harmonise_data(exposure_dat = exposure_data,
                               outcome_dat = outcome_data,
                               action = 3)
H_data <- filter(H_data, H_data$mr_keep == TRUE)

# Exclude H_data with Female-Male and Male-Female
H_data <- H_data %>%
  mutate(
    exposure_sex = str_extract(exposure, "(Female|Male)$"),
    outcome_sex  = str_extract(outcome, "(female|male)$") %>% str_to_title()
  )
H_data <- H_data %>%
  filter(exposure_sex == outcome_sex)

# Perform RadialMR analysis
H_data$exposure_outcome <- paste0(H_data$exposure, "_", H_data$outcome)
exposure_outcome_list <- unique(H_data$exposure_outcome)

data_all <- data.frame(matrix(ncol = ncol(H_data), nrow = 0))
colnames(data_all) <- colnames(H_data)

      for (i in 1:length(exposure_outcome_list)) {
        data <- filter(H_data, H_data$exposure_outcome == exposure_outcome_list[i])

        data_format <- format_radial(
          BXG = data$beta.exposure,
          BYG = data$beta.outcome,
          seBXG = data$se.exposure,
          seBYG = data$se.outcome,
          RSID = data$SNP
        )

        result <- ivw_radial(data_format, alpha = 0.05, weights = 3)

        if (!is.character(result$outliers) || result$outliers != "No significant outliers") {
          data <- anti_join(data, result$outliers, "SNP")
        }

        data_all <- rbind(data_all, data)
      }

H_data <- data_all

write.csv(H_data, file.path(output_dir, "H_data_for_MR_analysis.csv"), row.names = FALSE)

# Perform MR analysis (primary estimator: IVW)
mr_results <- generate_odds_ratios(mr(H_data, method_list = c("mr_ivw")))
write.csv(mr_results, file.path(output_dir, "MR_result_IVW_method.csv"), row.names = FALSE)

# Perform sensitivity tests
het <- mr_heterogeneity(H_data)
write.csv(het, file.path(output_dir, "sensitivity_cochrans_q.csv"), row.names = FALSE)

pleio <- mr_pleiotropy_test(H_data)
write.csv(pleio, file.path(output_dir, "sensitivity_mr_egger_intercept.csv"), row.names = FALSE)

# filter significant results for Sensitivity test
mr_results <- mr_results %>%
  group_by(outcome) %>%              # AD female / AD male stratified fdr
  mutate(pval_fdr = p.adjust(pval, method = fdr_method)) %>%
  ungroup()

sig_res <- filter(mr_results, pval_fdr < fdr_threshold)
H_data$exposure_outcome <- paste0(H_data$exposure,
                                  "_",
                                  H_data$outcome)
sig_res$exposure_outcome <- paste0(sig_res$exposure,
                                   "_",
                                   sig_res$outcome)
H_data_sig_res <- semi_join(H_data, sig_res, "exposure_outcome")

# significant MR results (all MR methods)
set.seed(2025) # Reproducibility
mr_results <- generate_odds_ratios(mr(H_data_sig_res, method_list = c(
  "mr_ivw",
  "mr_egger_regression",
  "mr_weighted_median",
  "mr_weighted_mode")))
write.csv(mr_results, file.path(output_dir, "MR_results_significant_all_methods.csv"), row.names = FALSE)

# leave-one-out analysis for significant results
dat_leave <- mr_leaveoneout(H_data_sig_res)
write.csv(dat_leave, file.path(output_dir, "sensitivity_leave_one_out.csv"), row.names = FALSE)

# MR-PRESSO
H_datapressoall <- filter(H_data, mr_keep == 'TRUE')
H_datapressoall$exposure_outcome <- paste0(H_datapressoall$exposure,
                                           "_",
                                           H_datapressoall$outcome)

testnum <- unique(H_datapressoall$exposure_outcome)

sum <- list()

for (i in seq_along(testnum)) {

  H_datapresso <- dplyr::filter(
    H_datapressoall,
    exposure_outcome == testnum[i]
  )

  SummaryStats <- data.frame(
    Y_effect  = H_datapresso$beta.outcome,
    E1_effect = H_datapresso$beta.exposure,
    Y_se      = H_datapresso$se.outcome,
    E1_se     = H_datapresso$se.exposure
  )

  a <- mr_presso(
    BetaOutcome = "Y_effect",
    BetaExposure = "E1_effect",
    SdOutcome = "Y_se",
    SdExposure = "E1_se",
    OUTLIERtest = TRUE,
    DISTORTIONtest = TRUE,
    data = SummaryStats,
    NbDistribution = mrpresso_nb_dist,
    SignifThreshold = mrpresso_signif
  )

  sum[[i]] <- data.frame(
    exposure  = H_datapresso$exposure[1],
    outcome   = H_datapresso$outcome[1],
    mainP     = a[["Main MR results"]][1, 6],
    correctP  = a[["Main MR results"]][2, 6],
    mainOR    = a[["Main MR results"]][1, 3],
    correctOR = a[["Main MR results"]][2, 3],
    mainsd    = a[["Main MR results"]][1, 4],
    correctsd = a[["Main MR results"]][2, 4],
    SNP       = length(a[["MR-PRESSO results"]][["Distortion Test"]][["Outliers Indices"]]),
    RSSobs    = a[["MR-PRESSO results"]][["Global Test"]][["RSSobs"]],
    globalP   = a[["MR-PRESSO results"]][["Global Test"]][["Pvalue"]]
  )
}

sum <- do.call(rbind, sum)

write.csv(sum, file.path(output_dir, "sensitivity_mrpresso.csv"), row.names = FALSE)
