# =============================================================================
# 04_motif_analysis / motifbreakR_ad_hormone.R
#
# Allele-specific motif disruption/creation for the prioritized AD risk SNPs in
# sex-hormone receptor motifs (AR, ESR1, ESR2, PGR).
#
# Two runs are performed with the same motif set:
#   1. score-based filtering:  threshold = 0.85,  filterp = FALSE
#   2. p-value filtering:      threshold = 1e-4,  filterp = TRUE, calculatePvalue()
#
# Required input columns: SNP (rsID), A1, A2 (effect / other allele).
# =============================================================================

# =========================
# USER CONFIGURATION
# =========================
input_file  <- "path/to/snps_with_alleles.csv"     # columns: SNP, A1, A2
output_dir  <- "path/to/output/motifbreakR/"

score_threshold <- 0.85          # threshold for the score-based run (filterp = FALSE)
p_threshold     <- 1e-4          # threshold for the p-value run (filterp = TRUE)
pct_threshold   <- 0.85          # keep pctRef >= pct_threshold or pctAlt >= pct_threshold
motif_gene_symbols <- c("AR", "ESR1", "ESR2", "PGR")
# =========================

library(motifbreakR)
library(BiocParallel)
library(SNPlocs.Hsapiens.dbSNP155.GRCh38)
library(BSgenome.Hsapiens.UCSC.hg38)
library(MotifDb)
library(Cairo)
library(data.table)

filter_mb_to_input_alleles <- function(mb_df, meta_df, a1_col, a0_col) {
  if (is.null(mb_df) || nrow(mb_df) == 0) {
    return(data.table())
  }

  mb_df[, REF := as.character(REF)]
  mb_df[, ALT := as.character(ALT)]

  keep_dt <- unique(meta_df[, .(
    SNP_id = SNP,
    A1 = get(a1_col),
    A0 = get(a0_col)
  )])

  out <- merge(as.data.table(mb_df), keep_dt, by = "SNP_id", all.x = TRUE)

  out <- out[
    (REF == A1 & ALT == A0) |
      (REF == A0 & ALT == A1)
  ]

  out[, c("A1", "A0") := NULL]
  out[]
}

granges_to_dt <- function(x) {
  if (length(x) == 0) {
    return(data.table())
  }
  as.data.table(as.data.frame(x, row.names = NULL))
}

df <- fread(input_file)

# keep only SNPs with an rsID
df <- df[grepl("^rs", SNP)]

# select top snps
snps <- df$SNP

variants <- snps.from.rsid(
  rsid = snps,
  dbSNP = SNPlocs.Hsapiens.dbSNP155.GRCh38,
  search.genome = BSgenome.Hsapiens.UCSC.hg38
)

motif_list <- lapply(motif_gene_symbols, function(gene_symbol) {
  subset(MotifDb, organism == "Hsapiens" & geneSymbol == gene_symbol)
})

motif_names <- unique(unlist(lapply(motif_list, names)))

hormone_motifs <- MotifDb[motif_names]

# -----------------------------
# 1) score-based results
# -----------------------------
mb_score_gr <- motifbreakR(
  snpList = variants,
  pwmList = hormone_motifs,
  threshold = score_threshold,
  filterp = FALSE,
  method = "ic",
  BPPARAM = SerialParam()
)

# select proper columns for output
if (length(mb_score_gr) != 0) {
  mb_score <- filter_mb_to_input_alleles(
    mb_df   = granges_to_dt(mb_score_gr),
    meta_df = df,
    a1_col  = "A1",
    a0_col  = "A2"
  )
  mb_score[, CHR := as.character(seqnames)]
  mb_score[, BP := start]
  mb_score[, motifPos_str := vapply(motifPos, function(x) paste(x, collapse = ","), character(1))]

  out <- mb_score[, .(
    SNP_id, CHR, BP,
    REF, ALT,
    geneSymbol,
    dataSource, providerName, providerId,
    seqMatch,
    motifPos_str,
    pctRef, pctAlt,
    scoreRef, scoreAlt,
    alleleDiff, alleleEffectSize,
    effect
  )]

  fwrite(out, file.path(output_dir, "ad_hormone_snps_mb_score.csv"))
}

# -----------------------------
# 2) p-value filter results
# -----------------------------
mb_P_gr <- motifbreakR(
  snpList = variants,
  pwmList = hormone_motifs,
  threshold = p_threshold,
  filterp = TRUE,
  method = "ic",
  BPPARAM = SerialParam()
)

if (length(mb_P_gr) != 0) {
  mb_P_gr <- calculatePvalue(mb_P_gr)
  mb_P <- filter_mb_to_input_alleles(
    mb_df   = granges_to_dt(mb_P_gr),
    meta_df = df,
    a1_col  = "A1",
    a0_col  = "A2"
  )

  keep_cols <- c(
    "SNP_id",
    "seqnames", "start",
    "REF", "ALT",
    "geneSymbol",
    "dataSource", "providerName", "providerId",
    "seqMatch", "motifPos",
    "pctRef", "pctAlt",
    "scoreRef", "scoreAlt",
    "Refpvalue", "Altpvalue",
    "alleleDiff", "alleleEffectSize",
    "effect", "pvalueEffect"
  )

  mb_P_out <- mb_P[, ..keep_cols]
  fwrite(mb_P_out, file.path(output_dir, "ad_hormone_snps_mb_P.csv"))
}

if (length(mb_score_gr) != 0 & length(mb_P_gr) != 0) {
   mb_P <- mb_P[pctRef >= pct_threshold | pctAlt >= pct_threshold]
   if(nrow(mb_P) != 0) {
     fwrite(mb_P, file.path(output_dir, "ad_hormone_snps_mb_P_filtered.csv"))
   }
}
