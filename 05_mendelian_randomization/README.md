# 05 - Mendelian randomization

## Purpose

Two-sample Mendelian randomization (MR) between sex hormone-related traits and Alzheimer's disease with the main estimator and a full set of sensitivity analyses, plus multivariable MR (MVMR) of SHBG accounting for BMI and the MRlap sample-overlap correction.

## Inputs

- **Instruments**: harmonised, LD-clumped instrument data for the sex hormone-related exposures, saved as an `.RData` object named `exposure_data` with columns `SNP`, `beta.exposure`, `se.exposure`, `effect_allele.exposure`, `other_allele.exposure`, `eaf.exposure`, `pval.exposure`, `samplesize.exposure`, `exposure`, `id.exposure`. Instruments were selected at P < 5e-8 and clumped (r2 = 0.001, 10,000 kb, EUR).
- **Outcome**: AD GWAS summary statistics saved as an `.RData` object named `trait_all` whose columns, in order, are chromosome, position, SNP, effect allele, other allele, effect allele frequency, beta, standard error, p-value, sample size and outcome label (they are renamed to the `*.outcome` columns by `mr_main_and_sensitivity.R`).
- **MVMR**: two exposure files with columns `exposure`, `SNP`, `beta.exposure`, `se.exposure`, `eaf.exposure`, `effect_allele.exposure`, `other_allele.exposure`, `pval.exposure`, `id.exposure` (BMI and SHBG), plus a PLINK binary reference panel (EUR) for clumping.
- **MRlap**: the MR IVW result table and the harmonised MR data table written by `mr_main_and_sensitivity.R`, the hormone trait and AD summary statistics with the columns listed in `mrlap.R`, and the GenomicSEM LD reference files (`eur_w_ld_chr/`) and HapMap3 SNP list (`w_hm3.snplist`), which are distributed separately by the GenomicSEM/MRlap documentation.

## Main scripts

| Script | Purpose |
|---|---|
| `mr_main_and_sensitivity.R` | Harmonises instruments and outcome, removes RadialMR outliers, runs IVW as the main analysis, then Cochran's Q, MR-Egger intercept, MR-Egger / weighted median / weighted mode for the FDR-significant results, leave-one-out, and MR-PRESSO |
| `mvmr.R` | Multivariable MR for SHBG and BMI on sex-stratified AD |
| `mrlap.R` | Sample-overlap correction with MRlap for the FDR-significant IVW results |

## Outputs

- `H_data_for_MR_analysis.csv` - harmonised data after RadialMR outlier removal
- `MR_result_IVW_method.csv` - IVW results (primary analysis)
- `sensitivity_cochrans_q.csv`, `sensitivity_mr_egger_intercept.csv` - heterogeneity and directional pleiotropy tests
- `MR_results_significant_all_methods.csv` - IVW, MR-Egger, weighted median and weighted mode for FDR-significant results
- `sensitivity_leave_one_out.csv` - leave-one-out estimates
- `sensitivity_mrpresso.csv` - MR-PRESSO global test, distortion test and outlier counts
- `mvdat.RData`, `mvmr_res.csv` - MVMR input data and results
- `mr_lap_res.xlsx` - MRlap corrected effect estimates
- one `<exposure>_<outcome>.Rdata` file per trait pair analysed by MRlap, holding the saved MRlap result object for that pair

## Key methodological parameters

- Instrument selection: P < 5e-8; clumping r2 = 0.001, window 10,000 kb, EUR reference panel; instruments restricted to the matched sex stratum (female exposure vs female outcome, male exposure vs male outcome).
- Instrument strength: F-statistic from `add_rsq` (`(N - 2) * rsq / (1 - rsq)`).
- RadialMR outlier removal: `format_radial` with beta/se pairs, `ivw_radial(alpha = 0.05, weights = 3)`.
- Main analysis: inverse-variance weighted (IVW) MR; odds ratios generated with `generate_odds_ratios`.
- Multiple testing: Benjamini-Hochberg FDR applied to the IVW p-values within each outcome.
- Sensitivity analyses (FDR-significant results): MR-Egger, weighted median, weighted mode, Cochran's Q (`mr_heterogeneity`), MR-Egger intercept (`mr_pleiotropy_test`), leave-one-out (`mr_leaveoneout`), MR-PRESSO (`OUTLIERtest = TRUE`, `DISTORTIONtest = TRUE`, `NbDistribution = 1000`, `SignifThreshold = 0.05`).
- Random seed: `set.seed(2025)` before the multi-method MR run and `set.seed(2025 + i)` per trait pair in the MRlap loop.
- MVMR: harmonisation with `harmonise_strictness = 3`, `mv_basic`, `mv_ivw`, `mv_multiple` and `mv_residual`, `pval_threshold = 5e-8`.
- MRlap: `do_pruning = FALSE` with `user_SNPsToKeep` set to the SNPs used in the corresponding MR analysis.

## Notes

- Plotting code (scatter plots, forest plots, leave-one-out plots) is not included; the leave-one-out estimates are written as a table instead.
- The expected column order of the outcome `.RData` object is documented above because the script renames the outcome columns positionally, as in the original analysis code.
- The GenomicSEM package required by MRlap must be installed separately (`devtools::install_github("GenomicSEM/GenomicSEM")`).
