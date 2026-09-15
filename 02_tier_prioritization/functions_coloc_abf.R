# =============================================================================
# Helper function for coloc_abf.R
#
# NOTE: the loop over the columns of the grouped locus table is kept exactly as
# in the analysis code (`for (j in 1:22)`). Verify that this range matches the
# column layout of your grouped locus table before re-running.
# =============================================================================

run_coloc_analysis <- function(sug_loci_results, AD_Meta_female, AD_Meta_male, AD_Meta_all, new_horm_data_list, PP4_results) {
  # Initialize PP4_results to store the results

  for (i in 1:nrow(sug_loci_results)) {
    snp <- sug_loci_results[i, ]
    print(paste("Processing SNP", i, "of", nrow(sug_loci_results), ":", snp$SNP))

    for (j in 1:22) {
      name <- names(sug_loci_results)[j]
      print(paste("Processing column:", name))
      if (!(snp[, j] < 0.05)) {
        print(paste("Skipping column:", name, "for SNP", snp$SNP))
        PP4_results[i, j] <- NA
        next
      }

      if (grepl("Meno", name)) {
        horm_data <- new_horm_data_list[["Meno"]]
      } else if (grepl("Mena", name)) {
        horm_data <- new_horm_data_list[["Mena"]]
      } else if (grepl("voice", name)) {
        horm_data <- new_horm_data_list[["voice"]]
      } else if (grepl("BTF", name)) {
        horm_data <- new_horm_data_list[["BTF"]]
      } else if (grepl("BTM", name)) {
        horm_data <- new_horm_data_list[["BTM"]]
      } else if (grepl("TTF", name)) {
        horm_data <- new_horm_data_list[["TTF"]]
      } else if (grepl("TTM", name)) {
        horm_data <- new_horm_data_list[["TTM"]]
      } else if (grepl("SBF", name)) {
        if (grepl("BMI", name)) {
          horm_data <- new_horm_data_list[["SBFBMI"]]
        } else {
          horm_data <- new_horm_data_list[["SBF"]]
        }
      } else if (grepl("SBM", name)) {
        if (grepl("BMI", name)) {
          horm_data <- new_horm_data_list[["SBMBMI"]]
        } else {
          horm_data <- new_horm_data_list[["SBM"]]
        }
      } else {
        horm_data <- NULL
      }

      if (grepl("female", name)) {
        AD_data <- AD_Meta_female
      } else if (grepl("male", name)) {
        AD_data <- AD_Meta_male
      } else if (grepl("all", name)) {
        AD_data <- AD_Meta_all
      } else {
        AD_data <- NULL
      }

      chr <- snp$CHR
      bp <- snp$BP
      range <- 1000000

      if(!snp$SNP %in% AD_data$SNP | !snp$SNP %in% horm_data$SNP) {
        print(paste("SNP", snp$SNP, "not found in AD or hormone data"))
        PP4 <- 0
        PP4_results[i, j] <- PP4
        next
      }

      # Filter AD data and hormone data within range
      AD_data_sub <- AD_data %>%
        dplyr::filter(CHR == chr & BP >= (bp - range) & BP <= (bp + range)) %>%
        mutate(varbeta = SE^2, id = paste(CHR, BP, A1, A2, sep = ":"))

      horm_data_sub <- horm_data %>%
        dplyr::filter(CHR == chr & BP >= (bp - range) & BP <= (bp + range)) %>%
        mutate(varbeta = SE^2, id12 = paste(CHR, BP, A1, A2, sep = ":"), id21 = paste(CHR, BP, A2, A1, sep = ":"))

      AD_data_sub <- AD_data_sub[AD_data_sub$BETA != 0,]
      horm_data_sub <- horm_data_sub[horm_data_sub$BETA != 0,]
      AD_data_sub <- AD_data_sub[AD_data_sub$varbeta != 0,]
      horm_data_sub <- horm_data_sub[horm_data_sub$varbeta != 0,]

      # Remove SNPs with no rsid
      AD_data_sub <- AD_data_sub[!is.na(AD_data_sub$SNP),]
      horm_data_sub <- horm_data_sub[!is.na(horm_data_sub$SNP),]
      #for duplicated SNPs, keep the one with the smallest p-value, only for AD data
      AD_data_sub <- AD_data_sub[order(AD_data_sub$P),]
      AD_data_sub <- AD_data_sub[!duplicated(AD_data_sub$SNP),]

      # Merge datasets by id12 and id21
      merged_AD_data_sub_12 <- merge(AD_data_sub[, c("CHR", "BP", "SNP", "EAF", "BETA", "SE", "varbeta", "N", "id")],
                                     horm_data_sub[, c("id12", "P")],
                                     by.x = "id", by.y = "id12") %>%
        dplyr::select(-P, -id)

      merged_horm_data_sub_12 <- merge(AD_data_sub[, c("P", "id")],
                                       horm_data_sub[, c("CHR", "BP", "SNP", "EAF", "BETA", "SE", "varbeta", "N", "id12")],
                                       by.x = "id", by.y = "id12") %>%
        dplyr::select(-P, -id)

      matched_SNP <- unique(merged_AD_data_sub_12$SNP)
      horm_data_sub <- horm_data_sub %>% dplyr::filter(!SNP %in% matched_SNP)

      merged_AD_data_sub_21 <- merge(AD_data_sub[, c("CHR", "BP", "SNP", "EAF", "BETA", "SE", "varbeta", "N", "id")],
                                     horm_data_sub[, c("id21", "P")],
                                     by.x = "id", by.y = "id21") %>%
        dplyr::select(-P, -id)

      merged_horm_data_sub_21 <- merge(AD_data_sub[, c("P", "id")],
                                       horm_data_sub[, c("CHR", "BP", "SNP", "EAF", "BETA", "SE", "varbeta", "N", "id21")],
                                       by.x = "id", by.y = "id21") %>%
        dplyr::select(-P, -id) %>%
        mutate(BETA = -BETA, EAF = 1 - EAF)

      merged_AD_data_sub <- rbind(merged_AD_data_sub_12, merged_AD_data_sub_21)
      merged_horm_data_sub <- rbind(merged_horm_data_sub_12, merged_horm_data_sub_21)

      # Filter out SNPs where EAF differs by more than 0.1
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

      # Prepare datasets for coloc analysis
      D1 <- list(
        beta = merged_AD_data_sub_filtered$BETA,
        varbeta = merged_AD_data_sub_filtered$varbeta,
        type = "cc",
        snp = merged_AD_data_sub_filtered$SNP,
        position = merged_AD_data_sub_filtered$BP,
        N = merged_AD_data_sub_filtered$N,
        MAF = merged_AD_data_sub_filtered$EAF
      )

      D2 <- list(
        beta = merged_horm_data_sub_filtered$BETA,
        varbeta = merged_horm_data_sub_filtered$varbeta,
        type = "quant",
        snp = merged_horm_data_sub_filtered$SNP,
        position = merged_horm_data_sub_filtered$BP,
        N = merged_horm_data_sub_filtered$N,
        MAF = merged_horm_data_sub_filtered$EAF
      )

      # Perform coloc analysis
      coloc_res <- coloc.abf(dataset1 = D1, dataset2 = D2)
      PP4 <- coloc_res$summary["PP.H4.abf"]
      PP4_results[i, j] <- PP4
      print("PP4_results updated")
    }
  }

  # Return the result matrix
  return(PP4_results)
}
