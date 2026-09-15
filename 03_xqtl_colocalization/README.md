# xQTL colocalization

## Purpose

Prioritize genes and regulatory features at Tier-1 pleiotropic loci by colocalizing Alzheimer's disease GWAS signals with molecular QTLs.

## QTL types represented

- eQTL
- sQTL
- pQTL
- mQTL
- caQTL
- haQTL

## Methods

- `coloc.abf` (`bulk_colocalization_abf.R` and `functions/*_abf.R`)
- `coloc.susie` (`bulk_colocalization_susie.R` and `functions/*_susie.R`)
- +/- 1 Mb analysis window around the locus
- strong colocalization defined as PP4 >= 0.70 (Tier-1 threshold)
- SuSiE analyses require an LD reference panel in PLINK binary format
- SuSiE settings: `runsusie` is called with `max_iter = 10000`, looping over `coverage` values from 0.95 down to 0.05 in steps of -0.1 and lowering the coverage until both datasets have at least one credible set

The xQTL colocalization functions were adapted from the Belloy Lab xQTL analysis workflow.

Molecular QTL datasets are **not** distributed in this repository and must be obtained from the original resources (see `data/README.md`). The scripts expect a QTL dataset that has already been formatted to a consistent variant/molecular-trait layout and harmonised to the same genome build as the AD GWAS. There is no automated download or dataset-formatting pipeline here.

## Directory

| Path | Content |
|---|---|
| `xqtl_coloc_abf.R` | Entry point for the coloc.abf workflow (USER CONFIGURATION block, then the bulk driver) |
| `xqtl_coloc_susie.R` | Entry point for the coloc.susie workflow (adds the PLINK/LD reference configuration) |
| `bulk_colocalization_abf.R` | Bulk driver: reads the formatted QTL dataset and the Tier-1 locus table, selects the AD stratum, harmonises AD and QTL variants, dispatches on `qtl_type`, combines results and writes one results table |
| `bulk_colocalization_susie.R` | Same workflow using the SuSiE implementations |
| `functions/eqtl_{abf,susie}.R` | eQTL colocalization |
| `functions/sqtl_{abf,susie}.R` | sQTL colocalization |
| `functions/pqtl_{abf,susie}.R` | pQTL colocalization |
| `functions/mqtl_{abf,susie}.R` | mQTL colocalization (CpG-level) |
| `functions/caqtl_{abf,susie}.R` | caQTL colocalization (peak-level, `Peak`) |
| `functions/haqtl_{abf,susie}.R` | haQTL colocalization (peak-level, `peak`) |

The bulk drivers dispatch eQTL, sQTL and pQTL, which are the QTL types harmonised by the shared driver code. The mQTL, caQTL and haQTL function files are retained because they are the QTL-specific colocalization implementations used after dataset-specific harmonisation (their molecular-trait key columns are `CpG`, `Peak` and `peak` respectively); harmonisation for those datasets was performed outside this driver and is therefore not represented here.

## Inputs

- `qtl_file`: formatted QTL dataset (one row per variant x molecular trait), with columns including `CHR`, `BP`, `id12`, `id21`, `gene_id`, `pvalue`, `beta`, `se`, `N`, `EAF`.
- `tier1_loci_file`: Tier-1 locus table with columns `chrom`, `bp_38_start`, `bp_38_end`, `bp_38_med`, `locus_index`, `stratum`, `discovery`. `stratum` selects the AD stratum (`Female`, `Male` or `NonStrat`) and `discovery == "pleiotropy"` selects the wider start/end window definition.
- AD GWAS summary statistics per stratum (case-control), with `posID`, `CHR`, `BP`, `SNP`, `ALLELE1`, `ALLELE0`, `A1FREQ`, `BETA`, `SE`, `P`, `N_incl`.
- GENCODE basic annotation (GTF) for gene boundaries.
- For SuSiE: `plink_bin`, `ld_reference_prefix` (PLINK binary prefix of the LD reference panel; the chromosome number is appended) and `temp_dir` for the per-locus SNP list and LD matrix files.

## Outputs

One results table per QTL dataset: `<study>_<tissue>_<qtl_type>_abf_hg38.csv` or `<study>_<tissue>_<qtl_type>_susie_hg38.csv`, with one row per locus/molecular trait and the columns `CHROM`, `BP`, `dset`, `locus`, `gene_name`, `molecular_trait_id`, `PP0`-`PP4`, `qtl_type`, `method`, `hit1`, `hit2`, `n_snps`, `discovery`.

## Notes on the original implementations (differences preserved)

- The AD/QTL effect-allele-frequency agreement filter differs between the two drivers: 0.1 in `bulk_colocalization_abf.R` and 0.5 in `bulk_colocalization_susie.R`. This difference is original and has been kept.
- The ABF implementations are not uniform: eQTL and pQTL build a single merged dataset, whereas sQTL, mQTL, caQTL and haQTL restrict the AD and QTL datasets separately, and the QTL-side filters differ (for example haQTL does not filter on `QTL_P > 0`). These QTL-specific choices were kept as they are.
- eQTL, sQTL and pQTL (`*_abf.R`) replace `NaN` posterior probabilities with 0; mQTL, caQTL and haQTL (`*_abf.R`) read the values directly from `coloc_results$summary` without that replacement.
- In the SuSiE implementations, PP0-PP4 are reported for the credible-set pair with the largest PP4; the sQTL implementation instead takes the maximum of each PP column separately (`max(..., na.rm = TRUE)`) and reports `hit1`/`hit2` from the first row. The `is.na(PP4)` guard and the SuSiE convergence/credible-set checks also differ between the QTL types (for example, eQTL and caQTL check the credible-set lengths, pQTL detects the "prior variance unreasonably large" error and skips that gene). All of these differences are original and were kept.
- The locus-compare plots generated by the original scripts are not part of this repository: the plotting blocks, the PDF creation, the `LC_dir` output logic and the plot-only column `lc_path` were removed. The PP4 >= 0.70 threshold is retained as the Tier-1 significance criterion and should be applied to the `PP4` column of the results table.
- The original SuSiE implementations contained a branch selecting an ancestry-specific LD panel for one project-specific stratum set. That branch was removed; the LD reference is now defined by `ld_reference_prefix` in the entry-point script.
- The eQTL SuSiE implementation included here is `autosome_eQTL_function_susie.R`. The original driver pointed at an alternative ("optimized") version of the same file, which is not part of this publication set.
- In the original ABF driver the eQTL function was called with positional arguments while the function signature had one parameter more than the other QTL types, so the arguments were shifted (the locus-compare output directory would have been written into the `dset` column). The public driver calls every QTL-type function with named arguments (`merged`, `all_suggestive_snps`, `name_dset`) so that the arguments match the signatures; no other part of the call was changed.
