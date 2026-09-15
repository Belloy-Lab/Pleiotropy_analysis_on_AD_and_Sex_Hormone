# Data

No data are distributed with this repository. This folder exists only to document where the input datasets come from and how to obtain them.

All datasets are obtained from their original providers. Controlled-access datasets (for example, individual-level or summary-level data derived from biobanks and consortia with data use agreements) remain subject to their original access requirements and cannot be redistributed here.

| Dataset category | Availability |
|---|---|
| Alzheimer disease GWAS | Publicly downloadable summary statistics or controlled-access release from the original study (Belloy et al., 2024), distributed with and without UK Biobank. UK Biobank-derived data must be accessed under the UK Biobank data access terms. |
| Sex hormone-related GWAS | Publicly downloadable summary statistics from the original studies / GWAS Catalog accessions (age at menarche, age at menopause, age at voice breaking, bioavailable testosterone, estradiol, sex hormone-binding globulin, total testosterone). UK Biobank-derived traits require UK Biobank access. |
| QTL datasets | Publicly downloadable summary statistics or controlled-access release from the original QTL consortia (for example eQTL Catalogue, eQTLGen, BrainMeta/MetaBrain and other tissue QTL resources). Some QTL resources require registration or an access request from the original provider. |
| LD/reference genotype data | The pleioFDR 1000 Genomes Phase 3 European LD reference files are distributed with the pleioFDR software. Genotype reference panels and unrelated-individual keep lists used for LD and SuSiE analyses must be obtained from their original providers, and the GenomicSEM / MRlap LD reference files and HapMap3 SNP list from the GenomicSEM documentation. |

Example input filenames used in the scripts are placeholders. Replace them with paths to files you obtained yourself, and make sure the column names and genome build expected by each script match your data (see the README in each analysis folder).
