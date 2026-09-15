# 02 - Tier-based locus prioritization and AD-hormone colocalization

## Purpose

Prioritization of the pleiotropic AD x sex hormone loci discovered with conjFDR into three tiers, and the component analyses used for that prioritization:

1. LD grouping of the conjFDR loci and annotation of the nearest/overlapping gene (`ld_grouping.R`, `annotate_genes.R`)
2. Colocalization between AD and the sex hormone-related trait with `coloc.abf` (`coloc_abf.R`)
3. Colocalization between AD and the sex hormone-related trait with `coloc.susie` (`coloc_susie.R`)
4. AD sex-heterogeneity P values and the female/male AD effect-size ratio (`sex_heterogeneity.R`)
5. Female/male AD effect-size comparison against the no-UKB replication (`effect_size_difference.R`)

## Tier criteria used in the manuscript

The final tier assignment reported in the manuscript was derived from the component analyses below and was **partly curated manually**. The scripts here provide those component analyses; they do not implement a single automatic tier-classification algorithm.

- **Tier 3**: conjFDR < 0.05
- **Tier 2**: Tier 3, plus either AD sex heterogeneity P < 0.05, or a >1.5-fold difference in the AD effect size between sexes
- **Tier 1**: Tier 2, plus AD-sex hormone-related trait colocalization PP4 >= 0.70

Which script computes which component:

| Component | Script | Output |
|---|---|---|
| AD sex-heterogeneity P < 0.05 (Tier 2) | `sex_heterogeneity.R` (first block) | `sex_het_p.tsv` |
| >1.5-fold difference in AD effect size between sexes (Tier 2) | `sex_heterogeneity.R` (second block) | `sex_het_p_1.5fold.tsv` (value `0.04`, i.e. below the 0.05 threshold, marks the ratio > 1.5 in either direction) |
| Female/male AD effect-size comparison and agreement with the no-UKB replication | `effect_size_difference.R` | `tier_check_with_delta*.csv`, `beta_comparison_r2.csv` |
| AD x hormone-related trait colocalization PP4 >= 0.70 (Tier 1) | `coloc_abf.R` (cross-checked with `coloc_susie.R`) | `sug_coloc_PP4_results.csv`, `*_susie_results.csv` |

## Inputs

- pleioFDR output for every AD x trait combination, in the layout expected by `ld_grouping.R` (see the header of that script): `results_<trait>_<AD stratum>_allexclude/result.mat` and `results_<trait>_<AD stratum>_allexclude/<ADstratumShort>_<trait>_conjfdr_0.05_loci.csv`, with `<AD stratum>` in `ADMeta_female`/`ADMeta_male`/`ADMeta_all`.
- pleioFDR LD reference table `9545380.ref` (distributed with pleioFDR).
- UCSC `liftOver` executable and the `hg19ToHg38.over.chain.gz` chain file.
- PLINK 1.9 and a genotype reference panel in PLINK binary format plus a keep list of unrelated individuals of matched ancestry, for LD estimation.
- AD GWAS summary statistics (female, male and non-stratified) and processed sex hormone-related trait summary statistics in hg38 (columns: `CHR`, `BP`, `SNP`/`rs37`, `A1`, `A2`, `EAF`, `BETA`, `SE`, `P`, `N`).
- A curated consensus table of known AD risk loci for the annotation step, with columns `CHR`, `GRCh38_POS`, `rsID`, `Locus_consensus`, `Study`.
- GENCODE basic annotation (GTF, protein-coding genes) in hg38.

## Main scripts

| Script | Purpose |
|---|---|
| `ld_grouping.R` | Collects all AD x trait conjFDR results, builds the per-locus FDR matrix, adds alleles, lifts hg19 -> hg38 (liftOver) and groups loci within 1 Mb into LD groups (PLINK r2) |
| `functions_ld_grouping.R` | Helper functions for `ld_grouping.R` (all paths passed as arguments) |
| `annotate_genes.R` | Annotates each LD group with the overlapping / nearest protein-coding gene and classifies it as a known AD locus (with/without LD to the reported variant) or a novel locus |
| `functions_annotation.R` | Helper functions for `annotate_genes.R` |
| `coloc_abf.R` | AD x sex hormone-related trait colocalization with `coloc.abf` for every locus/trait pair with conjFDR < 0.05 |
| `functions_coloc_abf.R` | Helper function for `coloc_abf.R` |
| `coloc_susie.R` | AD x sex hormone-related trait colocalization with `coloc.susie` for a single locus/trait pair (called once per locus, e.g. from a job array) |
| `sex_heterogeneity.R` | Annotates each index variant of the grouped locus table with the AD sex-heterogeneity P value, and flags variants with a >1.5-fold female/male AD effect-size ratio (Tier-2 components) |
| `effect_size_difference.R` | Computes female, male and sex-difference effect-size deltas between the AD GWAS and its no-UKB replication, joins the per-variant PP4_max and sex-heterogeneity P values, selects one variant per LD group, and reports the corresponding R2 values |

## Outputs

- `sug_fdrmat_results.csv` - per-locus conjFDR matrix across all AD x trait combinations
- `sug_fdrmat_results_withallele.csv`, `sug_fdrmat_results_hg38.csv`, `sug_fdrmat_results_hg38_grouped.csv` - intermediate and LD-grouped locus tables
- `sug_fdrmat_results_hg38_grouped_filtered.csv` - long-format table of locus/trait pairs with conjFDR < 0.05
- `sug_fdrmat_results_hg38_grouped_annotation.csv` - locus table with gene annotation
- `sug_coloc_PP4_results.csv` - colocalization PP4 per locus and trait (coloc.abf)
- one `*_susie_results.csv` per locus/trait pair (coloc.susie)
- `sex_het_p.tsv` - AD sex-heterogeneity P value of each index variant (written into the trait columns of the grouped locus table)
- `sex_het_p_1.5fold.tsv` - same table with the value `0.04` marking variants whose female/male AD effect-size ratio exceeds 1.5
- `tier_check_with_delta.csv`, `tier_check_with_delta_best_PP4_by_group.csv`, `tier_check_with_delta_best_sex_het_by_group.csv` - candidate locus table with effect-size deltas, `PP4_max` and sex-heterogeneity P values, plus the per-LD-group selections
- `beta_comparison_r2.csv` - squared correlation between the AD GWAS and its no-UKB replication for the female, male and sex-difference effect sizes

## Key methodological parameters

- LD grouping: loci in LD with r2 > 0.01 (PLINK `--r2 yes-really --matrix`, `--keep-allele-order`, `--allow-no-sex`) are merged; LD is only computed for SNP pairs on the same chromosome less than 1 Mb apart.
- AD-hormone colocalization window: +/- 1 Mb around the index variant.
- SNP harmonisation before colocalization: `varbeta = SE^2`; SNPs with `BETA == 0` or `varbeta == 0` removed; duplicate AD SNPs reduced to the lowest P value; strand matching through both allele orders (`id12`/`id21`), with beta and EAF flipped for the reverse order; SNPs with a difference in effect allele frequency > 0.1 between the two datasets removed.
- `coloc.abf`: default coloc priors; AD data as `type = "cc"`, hormone trait data as `type = "quant"`.
- `coloc.susie`: single effect per locus (`runsusie`, `max_iter = 10000`); the coverage is lowered in steps of 0.1 from 0.95 to 0.05 until both datasets have at least one credible set; reported PP4 is the maximum PP4 across credible-set pairs.
- Locus annotation: SNPs within 1 Mb of a curated known AD risk locus are compared with that locus by LD (r2 > 0.1 = known locus with LD); otherwise the SNP is annotated by position (overlapping protein-coding gene, or nearest gene by start position).
- Gene/annotation labels are propagated across all SNPs of an LD group (`process_groups`).
- Tier 2, sex heterogeneity: the AD sex-heterogeneity P value of the index variant is written into the trait columns of the grouped locus table; variants absent from the sex-heterogeneity GWAS are set to 1 (not significant), and a variant not found in the sex-stratified AD GWAS in the ratio step is also set to 1.
- Tier 2, effect-size difference: the female/male AD effect-size ratio is computed in both directions and a ratio > 1.5 is flagged with the value `0.04`.
- `effect_size_difference.R` summarises the colocalization results with `PP4_max = max(PP4_max)` per variant and keeps the variant with the largest `PP4_max` per LD group; the sex-heterogeneity P value is summarised per variant with `min(SEX_HET_P)` and one variant per group is kept with `slice_max(sex_het_p_max)`; R2 is the squared Pearson correlation.

## Notes and known limitations

- The tier labels are not produced by these scripts. Tier 2 and Tier 1 membership was curated from the component analyses above (see the tier criteria section) and should be documented in the manuscript supplementary tables.
- The `for (j in 1:22)` loop in `coloc_abf.R` matches the layout of the grouped locus table written by `ld_grouping.R`, in which the trait columns (11 traits x 2 AD strata = 22) come first and the identifier columns `CHR`, `BP`, `SNP` and `Group` come last. It is correct for the manuscript trait set, but the count is hard-coded: update `22` (and `n_trait_columns` in `sex_heterogeneity.R`) if the number of trait columns changes.
- `coloc.abf` is called with AD data of `type = "cc"` without an explicit case proportion `s`, so coloc's default handling applies.
- In `functions_annotation.R`, if the PLINK LD calculation fails for a SNP pair the code falls back to an r2 of 0.05 (below the 0.1 threshold, i.e. treated as no LD to the known locus). This fallback is preserved from the analysis code.
- The two Tier-2 component scripts originally read the AD sex-heterogeneity P values from different AD GWAS releases (one including UK Biobank, one without). Confirm which release is intended before re-running, and point both scripts at the same file.
- `effect_size_difference.R` keeps one variant per LD group with `slice_max(sex_het_p_max)`, i.e. the largest (least significant) sex-heterogeneity P value in the group, whereas the PP4 selection keeps the largest PP4. This is preserved from the original analysis script and should be confirmed.
