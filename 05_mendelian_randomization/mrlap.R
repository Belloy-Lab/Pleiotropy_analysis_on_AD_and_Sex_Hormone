# =============================================================================
# 05_mendelian_randomization / mrlap.R
#
# MRlap: sample-overlap correction of the MR estimates. MRlap is run for every
# IVW result that is significant after FDR correction, using the same SNPs that
# entered the corresponding two-sample MR analysis.
#
# MRlap requires the GenomicSEM LD reference files. Install GenomicSEM once:
#   devtools::install_github("GenomicSEM/GenomicSEM")
# =============================================================================

# =========================
# USER CONFIGURATION
# =========================
out_dir              <- "path/to/output/mrlap/"
mr_ivw_results_file  <- "path/to/output/mr/MR_result_IVW_method.csv"
h_data_file          <- "path/to/output/mr/H_data_for_MR_analysis.csv"   # after RadialMR exclusion

# Summary statistics in RData format:
#   hormone_sumstats_rdata -> object `hormone_trait` (columns CHR, BP, rsid, A1, A2, N, BETA, SE, P, trait)
#   ad_sumstats_rdata      -> object `trait_all`     (columns CHR, BP, SNP, ALLELE1, ALLELE0, N_incl, BETA, SE, P, trait)
hormone_sumstats_rdata <- "path/to/hormone_trait_id_matched.RData"
ad_sumstats_rdata      <- "path/to/ad_trait.RData"

# GenomicSEM reference files (obtained separately)
ld_dir    <- "path/to/GenomicSEMFiles/eur_w_ld_chr/"
hm3_file  <- "path/to/GenomicSEMFiles/w_hm3.snplist"

fdr_method <- "fdr"
# =========================

library(MRlap)
library(dplyr)
library(vroom)

# Generate analysis_list
mr_ivw_res <- read.csv(mr_ivw_results_file)
mr_ivw_res <- mr_ivw_res %>%
  group_by(outcome) %>%              # AD female / AD male stratified fdr
  mutate(pval_fdr = p.adjust(pval, method = fdr_method)) %>%
  ungroup()
mr_ivw_res <- filter(mr_ivw_res, mr_ivw_res$pval_fdr < 0.05)

analysis_list <- data.frame(hormone = mr_ivw_res$exposure,
                            AD = mr_ivw_res$outcome,
                            savename = paste0(mr_ivw_res$exposure, "_", mr_ivw_res$outcome, ".Rdata"))

H_data <- read.csv(h_data_file) #### generated from TwoSampleMR after RadialMR exclusion
H_data$exposure <- H_data$id.exposure
################################################################################
load(hormone_sumstats_rdata)
hormone_trait <- select(hormone_trait,
                        CHR,
                        BP,
                        rsid,
                        A1,
                        A2,
                        N,
                        BETA,
                        SE,
                        P,
                        trait)
colnames(hormone_trait) <- c("chr",
                     "pos",
                     "SNP",
                     "a1",
                     "a2",
                     "N",
                     "b",
                     "se",
                     "p",
                     "trait")
hormone_trait$a1 <- toupper(hormone_trait$a1)
hormone_trait$a2 <- toupper(hormone_trait$a2)
################################################################################


################################################################################
load(ad_sumstats_rdata)
AD_trait <- select(trait_all,
                   CHR,
                   BP,
                   SNP,
                   ALLELE1,
                   ALLELE0,
                   N_incl,
                   BETA,
                   SE,
                   P,
                   trait)
colnames(AD_trait) <- c("chr",
                       "pos",
                       "SNP",
                       "a1",
                       "a2",
                       "N",
                       "b",
                       "se",
                       "p",
                       "trait")
AD_trait$a1 <- toupper(AD_trait$a1)
AD_trait$a2 <- toupper(AD_trait$a2)
################################################################################

################################################################################
hormone_trait <- na.omit(hormone_trait)
AD_trait <- na.omit(AD_trait)

hormone_trait$chr <- as.numeric(hormone_trait$chr)
hormone_trait$pos <- as.numeric(hormone_trait$pos)
AD_trait$chr <- as.numeric(AD_trait$chr)
AD_trait$pos <- as.numeric(AD_trait$pos)
################################################################################

for (i in 1:nrow(analysis_list)){
  set.seed(2025 + i)

  save_name <- file.path(out_dir, analysis_list$savename[i])
  data1 <- filter(hormone_trait, trait == analysis_list$hormone[i])
  data2 <- filter(AD_trait, trait == analysis_list$AD[i])

  data1 <- distinct(data1, SNP, .keep_all = T)
  data2 <- distinct(data2, SNP, .keep_all = T)

  data2 <- semi_join(data2, data1, "SNP")
  data1 <- semi_join(data1, data2, "SNP")

  data1 <- arrange(data1, SNP)
  data2 <- arrange(data2, SNP)

  data1$trait <- NULL
  data2$trait <- NULL

  H_data_analysis <- filter(H_data, H_data$id.exposure == analysis_list$hormone[i] & H_data$id.outcome == analysis_list$AD[i])

  mrlap_res <- MRlap(exposure = data1,
                   exposure_name = paste0("exposure", i),
                   outcome = data2,
                   outcome_name = paste0("outcome", i),
                   ld = ld_dir,
                   hm3 = hm3_file,
                   do_pruning = FALSE,
                   user_SNPsToKeep = H_data_analysis$SNP)
  save(mrlap_res, file = save_name)
}

file_list <- dir(out_dir, pattern = "Rdata", full.names = TRUE)

mr_lap_res <- data.frame()
for (i in 1:length(file_list)) {
  load(file_list[i])
  MRcorrection <- mrlap_res[["MRcorrection"]]
  MRcorrection$IVs <- NULL
  df <- as.data.frame(MRcorrection)
  df$exposure = analysis_list$hormone[i]
  df$outcome = analysis_list$AD[i]
  mr_lap_res <- rbind(mr_lap_res, df)
}
# NOTE: the original analysis script had a stray trailing comma in this call
# (a syntax error); it has been removed so that the script parses. No analysis
# logic was changed.
mr_lap_res <- select(mr_lap_res,
                     exposure,
                     outcome,
                     m_IVs,
                     observed_effect,
                     observed_effect_se,
                     observed_effect_p,
                     corrected_effect,
                     corrected_effect_se,
                     corrected_effect_p)
writexl::write_xlsx(mr_lap_res, file.path(out_dir, "mr_lap_res.xlsx"))
