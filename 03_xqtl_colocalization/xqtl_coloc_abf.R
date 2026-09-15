# Publication version of the xQTL colocalization workflow
# supporting "Shared Risk Genes and Causal Relationships across
# Sex Hormone Related Traits and Alzheimer's Disease".
#
# Entry point for the bulk xQTL colocalization analysis with coloc.abf.
# Fill in the USER CONFIGURATION block below and run:
#
#   Rscript xqtl_coloc_abf.R
#
# The workflow is performed by bulk_colocalization_abf.R, which harmonises AD and QTL
# variants locus by locus and dispatches to the QTL-type specific implementations in
# functions/.

# =========================
# USER CONFIGURATION
# =========================
# Formatted QTL dataset (one row per variant x molecular trait) and the Tier-1 locus
# table that defines the loci/genes to test (columns: chrom, bp_38_start, bp_38_end,
# bp_38_med, locus_index, stratum, discovery).
qtl_file <- "path/to/formatted_qtl_dataset.txt.gz"
tier1_loci_file <- "path/to/tier1_loci_and_genes.csv"

# AD GWAS summary statistics per stratum
ad_female_file <- "path/to/ad/ad_female_gwas.txt.gz"
ad_male_file <- "path/to/ad/ad_male_gwas.txt.gz"
ad_nonstrat_file <- "path/to/ad/ad_nonstrat_gwas.txt.gz"

# GENCODE basic annotation (GTF) used to link molecular traits to genes, and output folder
gencode_file <- "path/to/reference/gencode.basic.annotation.gtf.gz"
output_dir <- "path/to/output/"

# Case proportion (cases / controls) of each AD stratum, passed to coloc as `s`
ad_female_s <- 59170/636266
ad_male_s <- 37435/507902
ad_nonstrat_s <- 95950/932843

# Labels of the QTL dataset; qtl_type also selects the colocalization implementation
study <- "QTL_study"
tissue <- "tissue"
qtl_type <- "eQTL"   # eQTL, sQTL or pQTL
# =========================

library(data.table)
library(dplyr)
library(rtracklayer)

# GENCODE annotation used to annotate QTL molecular traits with gene names/boundaries
ref_seq <- as.data.table(rtracklayer::import(gencode_file))

# Load one AD stratum and add the columns expected by the workflow:
#   Z  effect size / standard error
#   s  case proportion used by coloc for case-control data
load_ad_stratum <- function(file, s) {
  ad <- fread(file)
  ad$Z <- ad$BETA / ad$SE
  ad$s <- s
  return(ad)
}

Female_rb <- load_ad_stratum(ad_female_file, ad_female_s)
Male_rb <- load_ad_stratum(ad_male_file, ad_male_s)
NS_rb <- load_ad_stratum(ad_nonstrat_file, ad_nonstrat_s)

# Run the bulk colocalization workflow
source("bulk_colocalization_abf.R")
