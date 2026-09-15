# Shared Risk Genes and Causal Relationships across Sex Hormone Related Traits and Alzheimer's Disease

## Overview

This repository contains the analysis code supporting the manuscript "Shared Risk Genes and Causal Relationships across Sex Hormone Related Traits and Alzheimer's Disease". The code implements pleiotropy-informed conditional/conjunctional FDR discovery, tier-based locus prioritization with AD-hormone colocalization, xQTL colocalization, sex-hormone receptor motif analyses, and Mendelian randomization with sensitivity analyses. No data are distributed with this repository; all GWAS and QTL datasets must be obtained from their original sources (see [data/README.md](data/README.md)).

## Analysis workflow

1. **pleioFDR / conjFDR analysis** - genome-wide pleiotropy between Alzheimer's disease and sex hormone-related traits, conditional QQ plots and fold-enrichment: [01_pleiofdr](01_pleiofdr/README.md)
2. **Tier-based locus prioritization and AD-hormone colocalization** - LD grouping of conjFDR loci, gene annotation, and AD-hormone colocalization with `coloc.abf` and `coloc.susie`: [02_tier_prioritization](02_tier_prioritization/README.md)
3. **xQTL colocalization** - colocalization of AD risk loci with molecular QTL datasets across tissues: [03_xqtl_colocalization](03_xqtl_colocalization/README.md)
4. **Sex-hormone receptor motif analyses** - HOMER de novo motif enrichment and motifbreakR variant effect prediction for androgen, estrogen and progesterone receptor motifs: [04_motif_analysis](04_motif_analysis/README.md)
5. **Mendelian randomization** - two-sample MR with sensitivity analyses, multivariable MR, and sample-overlap correction: [05_mendelian_randomization](05_mendelian_randomization/README.md)

## Data availability

No data are included in this repository. All genome-wide association study (GWAS) and quantitative trait locus (QTL) summary statistics used in the manuscript must be obtained from their original sources, and controlled-access datasets remain subject to their original access requirements.

Data sources used by the code in this repository include:

- **Alzheimer's disease GWAS**: sex-stratified and non-stratified AD case-control GWAS summary statistics (Belloy et al., 2024), available with and without UK Biobank. UK Biobank-derived data are subject to UK Biobank access requirements.
- **Sex hormone-related GWAS**: age at menarche (23andMe-excluded release of Kentistou et al., 2024); age at menopause (reprogen, Ruth et al., 2021); age at voice breaking; bioavailable testosterone (GCST90012102, GCST90012103, GCST90012104); estradiol (GCST90012105); sex hormone-binding globulin with and without BMI adjustment (GCST90012106-GCST90012111); total testosterone (GCST90012112, GCST90012113, GCST90012114).
- **Negative-control traits**: hair colour GWAS (UK Biobank-derived), used to assess specificity of the pleiotropy signal.
- **Gene annotation and reference data**: GENCODE basic annotation (protein-coding genes); 1000 Genomes Phase 3 European LD reference (`9545380.ref`, distributed with pleioFDR) and a 9545380-SNP LD matrix reference distributed with pleioFDR; UCSC `liftOver` chain files for hg19/hg38 conversion; PLINK binary genotype reference panels and the accompanying unrelated-individual keep lists used for LD estimation.
- **QTL datasets**: bulk and other molecular QTL summary statistics (e.g. eQTL Catalogue, eQTLGen, BrainMeta/MetaBrain, and additional tissue QTL resources) as summarized in the manuscript, subject to the terms of the original data providers.

## Software

The analyses were performed with the following software. Specific versions should be recorded from the environment used for the final analysis; no version numbers are asserted here.

- **pleioFDR** (MATLAB) - conditional and conjunctional FDR; install the upstream implementation separately and add it to the MATLAB path, then run the wrapper in [01_pleiofdr](01_pleiofdr/README.md)
- **MATLAB** - runs the pleioFDR wrapper
- **R** packages: `data.table`, `dplyr`, `tidyr`, `stringr`, `reshape2`, `vroom`, `R.matlab`, `rtracklayer`, `Cairo`, `coloc`, `susieR`, `TwoSampleMR`, `MRPRESSO`, `RadialMR`, `ieugwasr`, `plinkbinr`, `MRlap` (with `GenomicSEM`), `writexl`
- **Bioconductor** packages (motif analysis): `motifbreakR`, `MotifDb`, `BSgenome.Hsapiens.UCSC.hg38`, `SNPlocs.Hsapiens.dbSNP155.GRCh38`, `BiocParallel`
- **PLINK** 1.9 - LD estimation (LD grouping, SuSiE LD matrices) and instrument clumping
- **liftOver** (UCSC) - hg19 to hg38 coordinate conversion
- **HOMER** (`findMotifs.pl`) - de novo transcription factor motif enrichment

## Citation

Manuscript in preparation. Authors, journal and DOI will be added here once the final bibliographic details are available.
