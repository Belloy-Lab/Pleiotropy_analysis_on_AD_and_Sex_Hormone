#!/usr/bin/env bash
#
# 04_motif_analysis / homer_find_motifs.sh
#
# De novo motif enrichment with HOMER (findMotifs.pl) for a foreground gene set
# against a background gene set, using a promoter window around the TSS and the
# HOCOMOCO v11 motif collection.
#
# The foreground and background gene lists (one gene symbol per line) must be
# prepared by the user; no gene lists are distributed with this repository.

set -euo pipefail

# =========================
# USER CONFIGURATION
# =========================
input_genes="path/to/foreground_genes.txt"        # foreground gene list (one gene symbol per line)
background_genes="path/to/background_genes.txt"   # background gene list (one gene symbol per line)
output_dir="path/to/output/homer_motifs/"         # output folder for findMotifs.pl
genome="human"                                    # HOMER genome assembly
motif_set="HOCOMOCOv11"                           # motif collection (-mset)
promoter_start=-1000                              # window start relative to the TSS (-start)
promoter_end=300                                  # window end relative to the TSS (-end)
score_threshold=30                                # motif score threshold (-S)
findmotifs_exe="findMotifs.pl"                    # path to findMotifs.pl if not on the PATH
# =========================

mkdir -p "$output_dir"

"$findmotifs_exe" \
  "$input_genes" \
  "$genome" \
  "$output_dir" \
  -bg "$background_genes" \
  -mset "$motif_set" \
  -start "$promoter_start" \
  -end "$promoter_end" \
  -S "$score_threshold"
