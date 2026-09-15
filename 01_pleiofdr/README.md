# 01 - pleioFDR / conjFDR analysis

## Purpose

Genome-wide pleiotropy analysis between Alzheimer's disease (AD) and sex hormone-related traits (and negative-control traits) with pleioFDR. The primary reported analysis is conjunctional FDR (conjFDR) at FDR < 0.05; conditional FDR, conditional QQ plots and fold-enrichment curves are produced by the same pleioFDR run.

## Inputs

- GWAS summary statistics for each trait, converted to pleioFDR `.mat` format. Use the upstream pleioFDR `sumstats.py` converter (`csv` -> `zscore` -> `mat`) with the pleioFDR reference file `9545380.ref` to build the `.mat` traits.
- A pleioFDR LD reference `.mat` file (for example `ref9545380_1kgPhase3eur_LDr2p1.mat`).
- The upstream pleioFDR MATLAB implementation (not redistributed here). Download it separately and point the `mlibrary` option in the config file to the folder containing `pleioOpt.m`, `TextConfig.m` and the pleioFDR analysis functions.

Traits analysed in the manuscript: AD (female, male, and non-stratified), age at menarche, age at menopause, age at voice breaking, bioavailable testosterone (female/male/all, GCST90012102-GCST90012104), estradiol (male, GCST90012105), sex hormone-binding globulin with and without BMI adjustment (female/male/all, GCST90012106-GCST90012111), total testosterone (female/male/all, GCST90012112-GCST90012114), and hair colour (female/male) as negative-control traits.

## Main scripts

- `run_pleiofdr.m` - wrapper: reads the config file, loads the LD reference, applies the exclusion regions and runs the upstream pleioFDR analysis.
- `config_pleiofdr_template.txt` - template config file. Copy it to `config_pleiofdr.txt` and edit the paths and trait files.

## Outputs

Written to the folder given by `outputdir` in the config file: pleioFDR result tables (`.mat` result objects, conjFDR/condFDR loci tables such as `<AD stratum>_<trait>_conjfdr_0.05_loci.csv`), conditional QQ plots, fold-enrichment plots and Manhattan plots.

## Key methodological parameters

Parameters below are the ones used in the manuscript config file (`config_pleiofdr_template.txt`); do not change them without a documented reason.

- `stattype=conjfdr` (recommended `fdrthresh` 0.05 for conjFDR, 0.01 for condFDR)
- `fdrthresh=0.05`, `pthresh=1` (no filtering on the Fisher combined statistic)
- `randprune=true`, `randprune_n=500`
- Exclusion regions (hg19): `[6 25119106 33854733; 8 7200000 12500000; 17 40000000 47000000]` (MHC, chromosome 8 inversion region, APOE region)
- `mafthresh=0.005`; ambiguous A/T and C/G SNPs excluded (`exclude_ambiguous_snps=true`)
- `exclude_from_discovery=true`
- Genomic control: `perform_gc=true`, `randprune_gc=true`, in-house (non-standard) genomic correction `use_standard_gc=false`
- Plot customisation: `manh_fontsize_genenames=12`, `manh_yspace=0.75`, `manh_ymargin=0.5`

## Notes

- The upstream pleioFDR MATLAB code (including the analysis/figure functions called by this wrapper) is not redistributed in this repository; see the pleioFDR project for its own licence and citation.
- The wrapper expects the pleioFDR analysis entry point to be on the MATLAB path; the config option `mlibrary` is added to the path automatically.
