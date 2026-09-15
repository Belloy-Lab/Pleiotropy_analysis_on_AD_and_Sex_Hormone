# 04 - Sex-hormone receptor motif analyses

## Purpose

Two complementary motif analyses of the prioritized AD risk variants:

1. **HOMER** - de novo motif enrichment in a promoter window around the prioritized genes, using a background gene set (`homer_find_motifs.sh`).
2. **motifbreakR** - prediction of allele-specific disruption/creation of sex-hormone receptor motifs (androgen receptor AR, estrogen receptors ESR1/ESR2, progesterone receptor PGR) for the prioritized SNPs (`motifbreakR_ad_hormone.R`).

## Inputs

### HOMER

- A foreground gene list (one gene symbol per line). The manuscript used a female-specific gene set.
- A background gene list (one gene symbol per line), for example all protein-coding genes, or a matched background set. Both foreground and background promoter/sequence sets must be prepared by the user; no gene lists are distributed with this repository.
- A local HOMER installation (`findMotifs.pl`).

### motifbreakR

- A table of prioritized SNPs with columns `SNP` (rsID, `rs...`), `A1` and `A2` (effect and other allele as reported for that SNP). Additional columns are ignored.
- Bioconductor annotation packages: `SNPlocs.Hsapiens.dbSNP155.GRCh38`, `BSgenome.Hsapiens.UCSC.hg38`, `MotifDb`.

## Main scripts

| Script | Purpose |
|---|---|
| `homer_find_motifs.sh` | Runs `findMotifs.pl` with a HOCOMOCO v11 motif set, a +/- promoter window and a score threshold, against a user-supplied background gene set |
| `motifbreakR_ad_hormone.R` | Runs `motifbreakR` twice for AR/ESR1/ESR2/PGR motifs: once with a score-based threshold (`filterp = FALSE`, threshold 0.85) and once with a p-value threshold (`filterp = TRUE`, threshold 1e-4) plus `calculatePvalue()`, and filters the results to the input alleles |

## Outputs

- HOMER: `homerResults/`, `knownResults/` and the usual `findMotifs.pl` output tables in the output folder.
- motifbreakR: `ad_hormone_snps_mb_score.csv` (score-filtered), `ad_hormone_snps_mb_P.csv` (p-value-filtered) and `ad_hormone_snps_mb_P_filtered.csv` (p-value-filtered, additionally restricted to variants with `pctRef >= 0.85` or `pctAlt >= 0.85`).

## Key methodological parameters

- HOMER: `-mset HOCOMOCOv11`, promoter window `-start -1000 -end 300` relative to the TSS, score threshold `-S 30`, background supplied with `-bg`.
- motifbreakR: `method = "ic"`; score-based run with `threshold = 0.85`, `filterp = FALSE`; p-value run with `threshold = 1e-4`, `filterp = TRUE` followed by `calculatePvalue()`; motif set restricted to human AR, ESR1, ESR2 and PGR motifs from `MotifDb`; results restricted to variants whose REF/ALT match the reported A1/A2 alleles in either orientation; final filtered table keeps variants with `pctRef >= 0.85` or `pctAlt >= 0.85`.
