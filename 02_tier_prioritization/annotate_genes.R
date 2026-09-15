# =============================================================================
# 02_tier_prioritization / annotate_genes.R
#
# Step 2 of the tier-based prioritization: annotate every LD group of conjFDR
# loci with
#   * the known AD risk locus it belongs to (if any), and whether the index
#     variant is in LD (r2 > 0.1) with the variant reported for that locus, or
#   * the overlapping protein-coding gene (GENCODE), or the nearest gene by
#     start position when the variant does not fall inside a gene.
# Annotation labels are then propagated to all SNPs of an LD group.
# =============================================================================

# =========================
# USER CONFIGURATION
# =========================
grouped_loci_file <- "path/to/output/sug_fdrmat_results_hg38_grouped.csv"
gtf_file          <- "path/to/reference/gencode.basic.annotation.gtf.gz"   # hg38, GENCODE
ad_risk_loci_file <- "path/to/reference/known_ad_risk_loci.txt"            # columns: CHR, GRCh38_POS, rsID, Locus_consensus, Study
plink_exe         <- "plink1.9"
ld_keep_file      <- "path/to/reference/unrelated_individuals.txt"
ld_bfile_prefix   <- "path/to/reference_genotypes/EUR_reference_panel_ch"  # chromosome number is appended
ld_work_dir       <- "path/to/work/ld_matrix/"
output_file       <- "path/to/output/sug_fdrmat_results_hg38_grouped_annotation.csv"
# =========================

library(data.table)
library(dplyr)
library(rtracklayer)

source("functions_annotation.R")

#------------------Load GTF File--------------------
gtf = rtracklayer::import(gtf_file)

#------------------Load Grouped FDR Results--------------------
sug_fdr_results <- fread(grouped_loci_file)

#------------------Process GTF File--------------------
ref_seq = as.data.frame(gtf)
ref_seq <- ref_seq[which(ref_seq$type=="gene" & ref_seq$gene_type=="protein_coding"),]

#------------------Load Known AD Risk Loci--------------------
AD_Risk_loci_consensus <- fread(ad_risk_loci_file)
sug_var_list <- sug_fdr_results %>% select(SNP,CHR,BP,Group)

#------------------Annotate Gene Names and Locus Information--------------------
sug_var_list$Gene_Name <- ""
sug_var_list$Gene_Name_near <- ""
sug_var_list$noval_locus <- ""
sug_var_list$known_locus_with_ld <- ""
sug_var_list$known_locus_no_ld <- ""
sug_var_list$Paper <- ""
sug_var_list$within_gene <- ""
sug_var_list$nearby_gene <- ""

#------------------Run Annotation Functions--------------------
updated_sug_var_list <- assign_gene_names(sug_var_list,
                                          AD_Risk_loci_consensus,
                                          ref_seq,
                                          plink_exe,
                                          ld_work_dir,
                                          ld_keep_file,
                                          ld_bfile_prefix)

processed_sug_var_list <- process_groups(updated_sug_var_list)

#------------------Save Processed Results--------------------
fwrite(processed_sug_var_list, output_file)
