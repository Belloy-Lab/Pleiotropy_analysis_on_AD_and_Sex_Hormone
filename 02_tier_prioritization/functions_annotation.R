# =============================================================================
# Helper functions for annotate_genes.R
#
# PLINK executable, genotype reference prefix, work directory and the keep list
# are passed in as arguments so that this file contains no machine-specific
# configuration.
# =============================================================================

assign_gene_names <- function(sug_var_list, AD_Risk_loci_consensus, ref_seq, plink_path, output_folder, keep_file, bfile_prefix) {

  # Start iterating over each unique group in the sug_var_list
  for (group in unique(sug_var_list$Group)) {

    print(paste("Processing group:", group))  # Print current group

    # Filter the SNPs belonging to the current group
    temp <- sug_var_list[sug_var_list$Group == group, ]
    print(paste("SNPs in this group:", paste(temp$SNP, collapse=", ")))  # Print the SNPs in the group

    # Extract the unique chromosome number for this group
    group_chr <- unique(temp$CHR)

    if (length(group_chr) != 1) {
      stop("Multiple chromosomes found in a single group, which is unexpected.")
    }

    # Check if group_chr exists in AD_Risk_loci_consensus
    if (!(group_chr %in% AD_Risk_loci_consensus$CHR)) {
      # If group_chr is not in the consensus data, set noval_locus to "Yes" and search for gene directly
      print(paste("Chromosome", group_chr, "not in AD_Risk_loci_consensus. Assigning noval_locus as 'Yes' and searching for gene."))

      sug_var_list$noval_locus[sug_var_list$Group == group] <- "Yes"

      # Filter the reference sequence data for genes on the same chromosome
      temp_ref_seq <- ref_seq[ref_seq$seqnames == paste0("chr", group_chr), ]

      # Iterate over SNPs in the group
      for (snp_index in 1:nrow(temp)) {
        snp_new <- temp$SNP[snp_index]
        bp_new <- temp$BP[snp_index]

        # Initialize a flag to track if a gene was found based on position
        gene_found <- FALSE

        # Search if SNP is within a gene
        for (i in 1:nrow(temp_ref_seq)) {
          if (temp_ref_seq$start[i] < bp_new & bp_new < temp_ref_seq$end[i]) {
            sug_var_list$Gene_Name_near[sug_var_list$SNP == snp_new] <- temp_ref_seq$gene_name[i]
            sug_var_list$within_gene[sug_var_list$SNP == snp_new] <- "Yes"
            gene_found <- TRUE
            print(paste("Gene found:", temp_ref_seq$gene_name[i]))
            break
          }
        }

        # If no gene was found within, assign the closest gene
        if (!gene_found) {
          closest_gene_index <- which.min(abs(temp_ref_seq$start - bp_new))
          closest_gene_name <- temp_ref_seq$gene_name[closest_gene_index]
          sug_var_list$Gene_Name_near[sug_var_list$SNP == snp_new] <- closest_gene_name
          sug_var_list$nearby_gene[sug_var_list$SNP == snp_new] <- "Yes"

          print(paste("Assigned closest gene:", closest_gene_name, "for SNP:", snp_new))
        }
      }

      # Continue to the next group since group_chr is not in the consensus
      next
    }

    # If the group contains only one SNP, execute this block
    if (nrow(temp) == 1) {
      print("Single SNP in group")

      # Filter the consensus data for SNPs on the same chromosome
      temp_concensus <- AD_Risk_loci_consensus %>% dplyr::filter(CHR == group_chr)

      # Extract the new SNP from the temp group
      snp_new <- temp$SNP
      bp_new <- temp$BP
      print(paste("Current SNP being processed:", snp_new))  # Print current SNP

      # Initialize a count variable to track number of LD matches
      count <- 0
      in_LD <- FALSE  # Initialize in_LD to track if LD was found

      # Check if the SNP is within 1Mb left or right of the known SNP
      bp_right <- bp_new + 1000000
      bp_left <- bp_new - 1000000
      in_region <- FALSE
      ld_result_pl <- NULL
      ld_result_ld <- NULL


      for (i in 1:nrow(temp_concensus)) {
        bp_known <- temp_concensus$GRCh38_POS[i]  # Get the known SNP from consensus data

        # Check if bp_known is not NA
        if (!is.na(bp_known)) {
          if (bp_left <= bp_known & bp_known <= bp_right) {
            # SNP is within 1Mb region of known SNP
            sug_var_list$noval_locus[sug_var_list$SNP == snp_new] <- "No"
            in_region <- TRUE
            in_file_chr <- paste0(bfile_prefix, group_chr)

            # Calculate LD (Linkage Disequilibrium)
            ld_result_pl <- tryCatch({
              calculate_r2_for_two_snps(snp_new, temp_concensus$rsID[i], plink_path, output_folder, keep_file, in_file_chr)
            }, error = function(e) {
              # If an error occurs, print a message and continue to next SNP
              message <- paste("Error with SNP:", snp_new, "and SNP:", temp_concensus$rsID[i], "-", e$message)
              print(message)  # Print error message for debugging
              return(NULL)  # Return NULL to safely continue the loop
            })

            if(is.null(ld_result_pl)) {
              # Calculate LD (Linkage Disequilibrium)
              ld_result_ld <- tryCatch({
                ld_result_ld = 0.05
              }, error = function(e) {
                # If an error occurs, print a message and return NULL
                print(paste("Error with SNP:", snp_new, "and SNP:", temp_concensus$rsID[i], "-", e$message))
                return(NULL)  # Return NULL to skip this SNP
              })
            }

            if(is.null(ld_result_pl) && is.null(ld_result_ld)) {
              print(paste("Skipping SNP:", snp_new, "due to null LD result."))
              next
            }

            if(!is.null(ld_result_pl)) {
              ld_result <- ld_result_pl
            } else {
              ld_result <- ld_result_ld
            }

            # Check if ld_result is valid
            if (!is.null(ld_result)) {
              if (!is.na(ld_result) && ld_result > 0.1) {
                # If LD (r2) is greater than 0.1, set known_locus_with_ld to "Yes"
                sug_var_list$known_locus_with_ld[sug_var_list$SNP == snp_new] <- "Yes"
                sug_var_list$Gene_Name[sug_var_list$SNP == snp_new] <- temp_concensus$Locus_consensus[i]
                sug_var_list$Paper[sug_var_list$SNP == snp_new] <- temp_concensus$Study[i]
                break  # Break out of the loop, no need to check further SNPs
              } else {
                # Otherwise, set known_locus_no_ld to "Yes"
                sug_var_list$known_locus_no_ld[sug_var_list$SNP == snp_new] <- "Yes"
                sug_var_list$Gene_Name[sug_var_list$SNP == snp_new] <- temp_concensus$Locus_consensus[i]
                sug_var_list$Paper[sug_var_list$SNP == snp_new] <- temp_concensus$Study[i]
              }
            } else {
              # Log when ld_result is NULL
              print(paste("Skipping SNP:", snp_new, "due to null LD result."))
            }
          }
        }
      }


      # If SNP is not within 1Mb of known SNP, set noval_locus to "Yes"
      if (!in_region) {
        sug_var_list$noval_locus[sug_var_list$SNP == snp_new] <- "Yes"
      }

      # If no LD matches were found, use genomic position (BP) to assign gene names
      if (!in_LD) {
        print("Searching by genomic position")

        # Filter the reference sequence data for genes on the same chromosome
        temp_ref_seq <- ref_seq[ref_seq$seqnames == paste0("chr", group_chr), ]

        # Initialize a flag to track if a gene was found based on position
        gene_found <- FALSE

        # Iterate through each gene in the reference sequence data
        for (i in 1:nrow(temp_ref_seq)) {
          # Check if the SNP is located between the start and end positions of a gene
          if (temp_ref_seq$start[i] < bp_new & bp_new < temp_ref_seq$end[i]) {
            sug_var_list$Gene_Name_near[sug_var_list$SNP == snp_new] <- temp_ref_seq$gene_name[i]
            sug_var_list$within_gene[sug_var_list$SNP == snp_new] <- "Yes"
            gene_found <- TRUE  # Set the flag to TRUE if a gene is found
            print(paste("Gene found:", temp_ref_seq$gene_name[i]))  # Print the gene found
            break  # Exit the loop if a gene is found
          }
        }

        # If no gene was found in the initial loop, assign the closest gene by start position
        if (!gene_found) {
          # Find the gene with the closest start position to the SNP
          closest_gene_index <- which.min(abs(temp_ref_seq$start - bp_new))
          closest_gene_name <- temp_ref_seq$gene_name[closest_gene_index]

          # Assign the closest gene to the SNP
          sug_var_list$Gene_Name_near[sug_var_list$SNP == snp_new] <- closest_gene_name
          sug_var_list$nearby_gene[sug_var_list$SNP == snp_new] <- "Yes"

          print(paste("Assigned closest gene:", closest_gene_name))  # Print closest gene assigned
        }
      }

      # Print the group, SNP, and count of LD matches
      print(paste("Final result for group:", group, "SNP:", snp_new, "LD matches:", count))
    }

    # If the group contains more than one SNP, execute this block
    if (nrow(temp) > 1) {
      print("Multiple SNPs in group")

      # Filter the consensus data for SNPs on the same chromosome
      temp_concensus <- AD_Risk_loci_consensus %>% filter(CHR == group_chr)

      # Iterate over each SNP in the current group
      for (snp_index in 1:nrow(temp)) {
        snp_new <- temp$SNP[snp_index]
        bp_new <- temp$BP[snp_index]
        print(paste("Processing SNP:", snp_new))  # Print current SNP

        # Initialize a count variable to track number of LD matches
        count <- 0
        in_LD <- FALSE  # Initialize in_LD to track if LD was found

        # Check if the SNP is within 1Mb left or right of any known SNP
        bp_right <- bp_new + 1000000
        bp_left <- bp_new - 1000000
        in_region <- FALSE
        ld_result_pl <- NULL
        ld_result_ld <- NULL

        for (i in 1:nrow(temp_concensus)) {
          bp_known <- temp_concensus$GRCh38_POS[i]  # Get the known SNP from consensus data

          # Check if bp_known is not NA
          if (!is.na(bp_known)) {
            if (bp_left <= bp_known & bp_known <= bp_right) {
              # SNP is within 1Mb region of known SNP
              sug_var_list$noval_locus[sug_var_list$SNP == snp_new] <- "No"
              in_region <- TRUE
              in_file_chr <- paste0(bfile_prefix, group_chr)

              # Calculate LD (Linkage Disequilibrium)
              ld_result_pl <- tryCatch({
                calculate_r2_for_two_snps(snp_new, temp_concensus$rsID[i], plink_path, output_folder, keep_file, in_file_chr)
              }, error = function(e) {
                # If an error occurs, print a message and continue to next SNP
                message <- paste("Error with SNP:", snp_new, "and SNP:", temp_concensus$rsID[i], "-", e$message)
                print(message)  # Print error message for debugging
                return(NULL)  # Return NULL to safely continue the loop
              })

              if(is.null(ld_result_pl)) {
                # Calculate LD (Linkage Disequilibrium)
                ld_result_ld <- tryCatch({
                  ld_result_ld = 0.05
                }, error = function(e) {
                  # If an error occurs, print a message and return NULL
                  print(paste("Error with SNP:", snp_new, "and SNP:", temp_concensus$rsID[i], "-", e$message))
                  return(NULL)  # Return NULL to skip this SNP
                })
              }

              if(is.null(ld_result_pl) && is.null(ld_result_ld)) {
                print(paste("Skipping SNP:", snp_new, "due to null LD result."))
                next
              }

              if(!is.null(ld_result_pl)) {
                ld_result <- ld_result_pl
              } else {
                ld_result <- ld_result_ld
              }

              if (!is.null(ld_result)) {
                # If LD (r2) is greater than 0.1, set known_locus_with_ld to "Yes"
                if (!is.na(ld_result) && ld_result > 0.1) {
                  sug_var_list$known_locus_with_ld[sug_var_list$SNP == snp_new] <- "Yes"
                  sug_var_list$Gene_Name[sug_var_list$SNP == snp_new] <- temp_concensus$Locus_consensus[i]
                  sug_var_list$Paper[sug_var_list$SNP == snp_new] <- temp_concensus$Study[i]
                  break  # Break out of the loop, no need to check further
                } else {
                  # Otherwise, set known_locus_no_ld to "Yes"
                  sug_var_list$known_locus_no_ld[sug_var_list$SNP == snp_new] <- "Yes"
                  sug_var_list$Gene_Name[sug_var_list$SNP == snp_new] <- temp_concensus$Locus_consensus[i]
                  sug_var_list$Paper[sug_var_list$SNP == snp_new] <- temp_concensus$Study[i]
                }
              }
            }
          }
        }

        # If SNP is not within 1Mb of known SNP, set noval_locus to "Yes"
        if (!in_region) {
          sug_var_list$noval_locus[sug_var_list$SNP == snp_new] <- "Yes"
        }

        # If no LD matches were found, use genomic position (BP) to assign gene names
        if (!in_LD) {
          print("Searching by genomic position")

          # Filter the reference sequence data for genes on the same chromosome
          temp_ref_seq <- ref_seq[ref_seq$seqnames == paste0("chr", group_chr), ]

          # Initialize a flag to track if a gene was found based on position.
          gene_found <- FALSE

          # Iterate through each gene in the reference sequence data
          for (i in 1:nrow(temp_ref_seq)) {
            # Check if the SNP is located between the start and end positions of a gene
            if (temp_ref_seq$start[i] < bp_new & bp_new < temp_ref_seq$end[i]) {
              sug_var_list$Gene_Name_near[sug_var_list$SNP == snp_new] <- temp_ref_seq$gene_name[i]
              sug_var_list$within_gene[sug_var_list$SNP == snp_new] <- "Yes"
              gene_found <- TRUE  # Set the flag to TRUE if a gene is found
              break  # Exit the loop if a gene is found
            }
          }

          # If no gene was found in the initial loop, assign the closest gene by start position
          if (!gene_found) {
            # Find the gene with the closest start position to the SNP
            closest_gene_index <- which.min(abs(temp_ref_seq$start - bp_new))
            closest_gene_name <- temp_ref_seq$gene_name[closest_gene_index]

            # Assign the closest gene to the SNP
            sug_var_list$Gene_Name_near[sug_var_list$SNP == snp_new] <- closest_gene_name
            sug_var_list$nearby_gene[sug_var_list$SNP == snp_new] <- "Yes"

            print(paste("Assigned closest gene:", closest_gene_name))  # Print closest gene assigned
          }
        }

        # Print the group, SNP, and count of LD matches
        print(paste("Final result for group:", group, "SNP:", snp_new, "LD matches:", count))
      }
    }
  }

  return(sug_var_list)
}

process_groups <- function(updated_sug_var_list) {

  # Ensure all "known_locus_with_ld" values are either "Yes" or "NA"
  updated_sug_var_list$known_locus_no_ld <- ifelse(
    updated_sug_var_list$known_locus_with_ld == "Yes",
    "",  # If 'known_locus_with_ld' is "Yes", set 'known_locus_no_ld' to an empty string
    updated_sug_var_list$known_locus_no_ld  # Else keep the original value
  )
  # Find the groups that contain more than one SNP (multiple rows)
  groups_to_process <- updated_sug_var_list %>%
    group_by(Group) %>%
    filter(n() > 1) %>%
    ungroup() %>%
    as.data.frame()  # Convert to data.frame to avoid Rle issues

  # Iterate over each unique group that has more than one row
  for (group in unique(groups_to_process$Group)) {

    # Filter the rows for the current group
    temp_group <- updated_sug_var_list %>% filter(Group == group) %>% as.data.frame()

    # Check if any row has `known_locus_with_ld` as "Yes"
    if (any(temp_group$known_locus_with_ld == "Yes")) {
      # Find the row with known_locus_with_ld == "Yes"
      row_with_ld <- temp_group[temp_group$known_locus_with_ld == "Yes", ][1, ]  # Take the first row if multiple

      # Update all rows in this group with the content from the row with known_locus_with_ld == "Yes"
      updated_sug_var_list[updated_sug_var_list$Group == group, c("Gene_Name", "Gene_Name_near","noval_locus", "known_locus_with_ld",
                                                                  "known_locus_no_ld", "Paper", "within_gene", "nearby_gene")] <-
        list(row_with_ld$Gene_Name, row_with_ld$Gene_Name_near ,row_with_ld$noval_locus, row_with_ld$known_locus_with_ld,
             row_with_ld$known_locus_no_ld, row_with_ld$Paper, row_with_ld$within_gene, row_with_ld$nearby_gene)

    } else {
      # Find the row with the BP closest to the median in the group
      median_bp <- median(temp_group$BP)
      row_closest_to_median <- temp_group[which.min(abs(temp_group$BP - median_bp)), ]

      # Update all rows in this group with the content from the row closest to the median BP
      updated_sug_var_list[updated_sug_var_list$Group == group, c("Gene_Name", "Gene_Name_near","noval_locus", "known_locus_with_ld",
                                                                  "known_locus_no_ld", "Paper", "within_gene", "nearby_gene")] <-
        list(row_closest_to_median$Gene_Name, row_closest_to_median$Gene_Name_near , row_closest_to_median$noval_locus, row_closest_to_median$known_locus_with_ld,
             row_closest_to_median$known_locus_no_ld, row_closest_to_median$Paper, row_closest_to_median$within_gene, row_closest_to_median$nearby_gene)
    }
  }

  return(updated_sug_var_list)
}

# r2 between two SNPs, estimated with PLINK from the genotype reference panel.
# `bfile` is the PLINK binary prefix for the relevant chromosome.
calculate_r2_for_two_snps <- function(snp1, snp2, plink_path, output_folder, keep_file, bfile) {

  # Create a temporary SNP list file named using snp1 and snp2
  common_snps_file <- file.path(output_folder, paste0("temp_", snp1, "_", snp2, ".txt"))

  # Write the two SNPs into the file
  write.table(data.frame(SNP = c(snp1, snp2)),
              file = common_snps_file,
              quote = FALSE, row.names = FALSE, col.names = FALSE)

  # Define output prefixes for PLINK commands to include snp1 and snp2
  out_prefix_data <- file.path(output_folder, paste0("temp_", snp1, "_", snp2, ".data"))
  out_prefix_ld <- file.path(output_folder, paste0("temp_", snp1, "_", snp2))

  # Generate a temporary PLINK binary dataset containing only the two SNPs
  plink_cmd1 <- paste(plink_path,
                      "--bfile", bfile,
                      "--allow-no-sex",
                      "--extract", common_snps_file,
                      "--keep-allele-order",
                      "--make-bed",
                      "--out", out_prefix_data)
  system(plink_cmd1)

  # Calculate LD (R2) matrix for the two SNPs
  plink_cmd2 <- paste(plink_path,
                      "--bfile", out_prefix_data,
                      "--allow-no-sex",
                      "--keep-allele-order",
                      "--keep", keep_file,
                      "--r2 yes-really --matrix",
                      "--out", out_prefix_ld)
  system(plink_cmd2)

  # Specify file paths for LD and BIM files
  ld_file <- paste0(out_prefix_ld, ".ld")
  bim_file <- paste0(out_prefix_data, ".bim")

  # Read BIM file to get SNP names
  bim <- fread(bim_file)
  snps_name <- bim$V2

  # Read LD matrix and convert to matrix format
  LD <- fread(ld_file)
  LD <- as.matrix(LD)
  colnames(LD) <- snps_name
  rownames(LD) <- snps_name

  # Check if both SNPs are present in the LD matrix
  if(!(snp1 %in% snps_name && snp2 %in% snps_name)) {
    stop("Error: SNPs not found in the LD file.")
    return(NULL)
  }

  # Extract R2 value for the given SNP pair
  r2_value <- LD[snp1, snp2]

  # Remove temporary files
  file.remove(common_snps_file)
  file.remove(paste0(out_prefix_data, ".bed"))
  file.remove(paste0(out_prefix_data, ".bim"))
  file.remove(paste0(out_prefix_data, ".fam"))
  file.remove(paste0(out_prefix_ld, ".nosex"))
  file.remove(paste0(out_prefix_ld, ".data.nosex"))
  file.remove(paste0(out_prefix_ld, ".ld"))

  return(r2_value)
}
