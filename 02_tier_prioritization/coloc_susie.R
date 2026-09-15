# =============================================================================
# 02_tier_prioritization / coloc_susie.R
#
# AD x sex hormone-related trait colocalization with coloc.susie for a single
# locus. The script is called once per locus/trait pair, for example from a job
# array:
#
#   Rscript coloc_susie.R --args <rsID> <CHR> <BP> <col_name>
#
# where <col_name> is the column name of the grouped locus table for that
# trait/AD stratum combination (e.g. "Meno_ADMeta_female_mat").
#
# Required columns
#   AD GWAS       : CHR, BP, SNP, ALLELE1, ALLELE0, A1FREQ, N_incl, BETA, SE, P
#   hormone trait : CHR, BP, A1, A2, EAF, BETA, SE, P, N and an rsid column
#                   (rs37; not used for age at voice breaking)
# =============================================================================

# =========================
# USER CONFIGURATION
# =========================
ad_female_file <- "path/to/ad/ad_female_gwas.txt.gz"
ad_male_file   <- "path/to/ad/ad_male_gwas.txt.gz"
ad_all_file    <- "path/to/ad/ad_all_gwas.txt.gz"

hormone_files <- list(
  Meno   = "path/to/traits/age_at_menopause_hg38.txt",
  Mena   = "path/to/traits/age_at_menarche_hg38.txt",
  voice  = "path/to/traits/age_at_voice_breaking_hg38.txt",
  BTF    = "path/to/traits/bioavailable_testosterone_female_hg38.txt",
  BTM    = "path/to/traits/bioavailable_testosterone_male_hg38.txt",
  TTF    = "path/to/traits/total_testosterone_female_hg38.txt",
  TTM    = "path/to/traits/total_testosterone_male_hg38.txt",
  SBFBMI = "path/to/traits/shbg_female_bmi_adjusted_hg38.txt",
  SBMBMI = "path/to/traits/shbg_male_bmi_adjusted_hg38.txt",
  SBF    = "path/to/traits/shbg_female_hg38.txt",
  SBM    = "path/to/traits/shbg_male_hg38.txt"
)

plink_exe       <- "plink1.9"
ld_bfile_prefix <- "path/to/reference_genotypes/EUR_reference_panel_ch"  # chromosome number is appended
ld_keep_file    <- "path/to/reference/unrelated_individuals.txt"
ld_work_dir     <- "path/to/work/coloc_susie/"                  # SNP lists and PLINK LD matrices
output_dir      <- "path/to/output/coloc_susie/"

locus_window   <- 1000000    # +/- window around the index variant
snp_list_prefix <- "coloc_snp"
ld_prefix       <- "coloc_snp_plink_ld_matrix"
# =========================

#------------------Load Libraries--------------------
library(data.table)
library(dplyr)
library(coloc)

#------------------Read Command Line Arguments--------------------
args <- commandArgs(trailingOnly = TRUE)

#------------------Assign Arguments to Variables--------------------
rs37 <- args[2]
CHR <- as.numeric(args[3])
BP <- as.numeric(args[4])
col_name <- args[5]

#------------------Print Information--------------------
cat("Running coloc analysis for:\n")
cat("rs37:", rs37, "\nCHR:", CHR, "\nBP:", BP, "\ncol_name:", col_name, "\n")

print("Starting loading AD data")
if (grepl("female", col_name)) {
    AD_data_path <- ad_female_file
} else if (grepl("male", col_name)) {
    AD_data_path <- ad_male_file
} else if (grepl("all", col_name)) {
    AD_data_path <- ad_all_file
}
print(AD_data_path)

print("Starting loading Horm data")
# Determine Horm_data_path based on col_name
if (grepl("Meno", col_name)){
    Horm_data_path <- hormone_files$Meno
} else if (grepl("Mena", col_name)){
    Horm_data_path <- hormone_files$Mena
} else if (grepl("voice", col_name)){
    Horm_data_path <- hormone_files$voice
} else if (grepl("BTF", col_name)){
    Horm_data_path <- hormone_files$BTF
} else if (grepl("BTM", col_name)){
    Horm_data_path <- hormone_files$BTM
} else if (grepl("TTF", col_name)){
    Horm_data_path <- hormone_files$TTF
} else if (grepl("TTM", col_name)){
    Horm_data_path <- hormone_files$TTM
}
if (grepl("SBF", col_name)){
    if (grepl("BMI", col_name)){
        Horm_data_path <- hormone_files$SBFBMI
    } else {
        Horm_data_path <- hormone_files$SBF
    }
} else if (grepl("SBM", col_name)){
    if (grepl("BMI", col_name)){
        Horm_data_path <- hormone_files$SBMBMI
    } else {
        Horm_data_path <- hormone_files$SBM
    }
}
print(Horm_data_path)

print("Starting reading AD and Horm data")

#------------------Load and Process Data--------------------
snp <- rs37
ad_data = fread(AD_data_path)
setnames(ad_data, old = c("ALLELE1", "ALLELE0", "A1FREQ", "N_incl","SNP"),
                new = c("A1", "A2", "EAF", "N", "rs38"))
ad_data$SNP <- ad_data$rs38
horm_data = fread(Horm_data_path)
if (!(grepl("voice", col_name))){
    horm_data$SNP <- horm_data$rs37
}
in_file_chr <- paste0(ld_bfile_prefix, CHR)
chr <- CHR
bp <- BP
range <- locus_window

#------------------Subset Data Around Locus--------------------
AD_data_sub <- ad_data %>%
    dplyr::filter(CHR == chr & BP >= (bp - range) & BP <= (bp + range)) %>%
    mutate(varbeta = SE^2) %>%
    mutate(id = paste(CHR, BP, A1, A2, sep = ":"))
horm_data_sub <- horm_data %>%
    dplyr::filter(CHR == chr & BP >= (bp - range) & BP <= (bp + range)) %>%
    mutate(varbeta = SE^2) %>%
    mutate(id12 = paste(CHR, BP, A1, A2, sep = ":"), id21 = paste(CHR, BP, A2, A1, sep = ":"))

AD_data_sub <- AD_data_sub[AD_data_sub$BETA != 0,]
horm_data_sub <- horm_data_sub[horm_data_sub$BETA != 0,]
AD_data_sub <- AD_data_sub[AD_data_sub$varbeta != 0,]
horm_data_sub <- horm_data_sub[horm_data_sub$varbeta != 0,]

#------------------Remove NA and Duplicated SNPs--------------------
AD_data_sub <- AD_data_sub[!is.na(AD_data_sub$SNP),]
horm_data_sub <- horm_data_sub[!is.na(horm_data_sub$SNP),]
#for duplicated SNPs, keep the one with the smallest p-value, only for AD data
AD_data_sub <- AD_data_sub[order(AD_data_sub$P),]
AD_data_sub <- AD_data_sub[!duplicated(AD_data_sub$SNP),]

#-------------------Merge Data for Coloc Analysis--------------------
merged_AD_data_sub_12 <- merge(AD_data_sub[, c("CHR", "BP", "SNP","EAF","BETA","SE","varbeta","N","id")], horm_data_sub[, c("id12", "P")], by.x = "id", by.y = "id12") %>% dplyr::select(-id, -P)
merged_horm_data_sub_12 <- merge(AD_data_sub[,c("P","id")], horm_data_sub[,c("CHR", "BP", "SNP","EAF","BETA","SE","varbeta","N","id12")], by.x = "id", by.y = "id12") %>% dplyr::select(-id, -P)

merged_AD_data_sub_21 <- merge(AD_data_sub[, c("CHR", "BP", "SNP","EAF","BETA","SE","varbeta","N","id")], horm_data_sub[, c("id21", "P")], by.x = "id", by.y = "id21") %>% dplyr::select(-id, -P)
merged_horm_data_sub_21 <- merge(AD_data_sub[,c("P","id")], horm_data_sub[,c("CHR", "BP", "SNP","EAF","BETA","SE","varbeta","N","id21")], by.x = "id", by.y = "id21") %>% dplyr::select(-id, -P) %>% mutate(BETA = -BETA, EAF = 1 - EAF)

merged_AD_data_sub <- rbind(merged_AD_data_sub_12, merged_AD_data_sub_21)
merged_horm_data_sub <- rbind(merged_horm_data_sub_12, merged_horm_data_sub_21)

#-------------------Filter SNPs by EAF Difference--------------------
combined_data <- full_join(
    merged_AD_data_sub %>% dplyr::select(SNP, EAF) %>% dplyr::rename(EAF_AD = EAF),
    merged_horm_data_sub %>% dplyr::select(SNP, EAF) %>% dplyr::rename(EAF_horm = EAF),
    by = "SNP"
)

filtered_snps <- combined_data %>%
    mutate(EAF_diff = abs(EAF_AD - EAF_horm)) %>%
    dplyr::filter(EAF_diff <= 0.1) %>%
    dplyr::select(SNP)

merged_AD_data_sub_filtered <- merged_AD_data_sub %>%
    dplyr::filter(SNP %in% filtered_snps$SNP)

merged_horm_data_sub_filtered <- merged_horm_data_sub %>%
    dplyr::filter(SNP %in% filtered_snps$SNP)

SNPs <- data.frame(merged_AD_data_sub_filtered[which(!duplicated(merged_AD_data_sub_filtered$SNP)), "SNP"])

n_of_snps <- nrow(SNPs)

#------------------Calculate LD Matrix Using PLINK--------------------
snp_list_file <- file.path(ld_work_dir, paste0(snp_list_prefix, rs37, col_name, ".txt"))
ld_out_prefix <- file.path(ld_work_dir, paste0(ld_prefix, rs37, col_name))

fwrite(SNPs, snp_list_file, sep = "\t", col.names = FALSE)
Sys.sleep(5)

command <- paste(plink_exe,
                   "--bfile", in_file_chr,
                   "--allow-no-sex",
                   "--extract", snp_list_file,
                   "--keep-allele-order",
                   "--keep", ld_keep_file,
                   "--make-bed --out", paste0(ld_out_prefix, ".data"),
                   sep= " ")
system(command)
Sys.sleep(5)

command <- paste(plink_exe,
                   "--bfile", paste0(ld_out_prefix, ".data"),
                   "--allow-no-sex",
                   "--keep-allele-order",
                   "--r --matrix",
                   "--out", ld_out_prefix,
                   sep= " ")
system(command)
Sys.sleep(5)

bim <- fread(paste0(ld_out_prefix, ".data.bim"))
snps <- unlist(bim$V2)

LD <- fread(paste0(ld_out_prefix, ".ld"))
LD <- as.matrix(LD)
colnames(LD) <- snps
rownames(LD) <- snps

#------------------Create SuSiE Input Lists--------------------
ADl <- as.list(merged_AD_data_sub_filtered[,c("BETA","varbeta","SNP","BP","EAF")])
names(ADl)[1:5] <- c("beta","varbeta","snp","position","MAF")
names(ADl$beta) <- ADl$snp
names(ADl$varbeta) <- ADl$snp
ADl$N <- median(merged_AD_data_sub_filtered$N, na.rm = TRUE)

matching_snps_ADl <- ADl$snp %in% snps
ADl$beta <- ADl$beta[matching_snps_ADl]
ADl$varbeta <- ADl$varbeta[matching_snps_ADl]
ADl$snp <- ADl$snp[matching_snps_ADl]
ADl$position <- ADl$position[matching_snps_ADl]
ADl$MAF <- ADl$MAF[matching_snps_ADl]
ADl$LD <- LD
ADl$type <- "cc"
cc <- ADl

common_snps <- intersect(cc$snp, colnames(cc$LD))

cc_index <- match(common_snps, cc$snp)
ld_index <- match(common_snps, colnames(cc$LD))

cc$beta <- cc$beta[cc_index]
cc$varbeta <- cc$varbeta[cc_index]
cc$snp <- cc$snp[cc_index]
cc$position <- cc$position[cc_index]
cc$MAF <- cc$MAF[cc_index]
cc$LD <- cc$LD[ld_index, ld_index]

#-------------------Create SuSiE Input Lists for Horm Data--------------------
HormL <- as.list(merged_horm_data_sub_filtered[,c("BETA","varbeta","SNP","BP","EAF")])
names(HormL)[1:5] <- c("beta","varbeta","snp","position","MAF")
names(HormL$beta) <- HormL$snp
names(HormL$varbeta) <- HormL$snp
HormL$N <- median(merged_horm_data_sub_filtered$N, na.rm = TRUE)

#------------------Match SNPs with LD Matrix--------------------
matching_snps_HormL <- HormL$snp %in% snps
HormL$beta <- HormL$beta[matching_snps_HormL]
HormL$varbeta <- HormL$varbeta[matching_snps_HormL]
HormL$snp <- HormL$snp[matching_snps_HormL]
HormL$position <- HormL$position[matching_snps_HormL]
HormL$MAF <- HormL$MAF[matching_snps_HormL]
HormL$LD <- LD
HormL$type <- "quant"
eq <- HormL

common_snps <- intersect(eq$snp, colnames(eq$LD))

eq_index <- match(common_snps, eq$snp)
ld_index <- match(common_snps, colnames(eq$LD))

eq$beta <- eq$beta[eq_index]
eq$varbeta <- eq$varbeta[eq_index]
eq$snp <- eq$snp[eq_index]
eq$position <- eq$position[eq_index]
eq$MAF <- eq$MAF[eq_index]
eq$LD <- eq$LD[ld_index, ld_index]

#------------------Run SuSiE and Coloc Analysis--------------------
coverage_values <- seq(0.95, 0.05, by = -0.1)
cs_length_S1 <- 0
cs_length_S2 <- 0
for (coverage in coverage_values) {
    S1 <- runsusie(cc, coverage = coverage, max_iter = 10000)
    cs_length_S1 <- length(S1$sets$cs)
    S2 <- runsusie(eq, coverage = coverage, max_iter = 10000)
    cs_length_S2 <- length(S2$sets$cs)
    if (cs_length_S1 > 0 & cs_length_S2 > 0) {
        break
    }
}

#-------------------Perform Coloc Analysis--------------------
out <- coloc.susie(S1, S2)

if (is.null(out$summary$PP.H4.abf)) {
    print("No coloc results")
    coloc_PP4 = 0
    coloc_hit1 = 0
    coloc_hit2 = 0
}else{
    print("Coloc results found")
    print(out$summary)
    coloc_PP4 = max(out$summary$PP.H4.abf)
    coloc_hit1 = out$summary$hit1[1]
    coloc_hit2 = out$summary$hit2[1]
}

#------------------Save Results--------------------
results_susie <- data.frame(chr = CHR, bp = BP, SNP = snp, pairs = col_name, hit1 = coloc_hit1, hit2 = coloc_hit2, PP4 = coloc_PP4, n_snps = n_of_snps, coverage = coverage, stringsAsFactors = FALSE)

file_name <- file.path(output_dir, paste0(rs37, "_", col_name, "_susie_results.csv"))
fwrite(results_susie, file_name, append = TRUE, col.names = TRUE)

#------------------Clean Up Temporary Files--------------------
suppressWarnings(file.remove(paste0(ld_out_prefix, ".data.bed")))
suppressWarnings(file.remove(paste0(ld_out_prefix, ".data.bim")))
suppressWarnings(file.remove(paste0(ld_out_prefix, ".data.fam")))
suppressWarnings(file.remove(paste0(ld_out_prefix, ".data.log")))
suppressWarnings(file.remove(paste0(ld_out_prefix, ".data.nosex")))
suppressWarnings(file.remove(paste0(ld_out_prefix, ".ld")))
suppressWarnings(file.remove(paste0(ld_out_prefix, ".log")))
suppressWarnings(file.remove(paste0(ld_out_prefix, ".nosex")))
suppressWarnings(file.remove(snp_list_file))
