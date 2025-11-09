### Load Required Libraries ###
# This script assumes all required packages have been installed.
# Use install.packages("package_name") or BiocManager::install("package_name") if any are missing.

library(phyloseq)
library(microbiomeMarker) # For LEfSe analysis
library(dplyr)
library(data.table)       # For the fread function
library(tidyr)
library(ape)


#' Run LEfSe Analysis
#'
#' This function encapsulates the entire LEfSe analysis workflow, from reading
#' the input data to performing the analysis.
#'
#' @param input_path Path to the input CSV file.
#' @param lda_cutoff LDA score cutoff for LEfSe analysis. Default is 3.
#' @param strict Strictness of the LEfSe analysis. Default is "1".
#' @param transform The transformation to apply to the data. Default is "log10p".
#' @param abundance_filter_threshold Threshold for filtering taxa by abundance. Default is 1000.
#' @param sample_min_threshold Minimum number of samples a taxon must be present in. Default is 3.
#' @param curv Whether to use a curved LDA model. Default is TRUE.
#' @return A `microbiomeMarker` object containing the LEfSe results.
run_lefse_analysis <- function(input_path = "data/processed/ALL-update.csv",
                               lda_cutoff = 3,
                               strict = "1",
                               transform = "log10p",
                               rarefy_seed = 394582,
                               abundance_filter_threshold = 1000,
                               sample_min_threshold = 3,
                               multigrp_strat = FALSE,
                               curv = FALSE) {

  # 1. Read the data using data.table for speed
  raw_data_in <- fread(input_path)

  # 2. Separate sample metadata
  sample_info_df <- raw_data_in %>%
    select(index, Group) %>%
    as.data.frame()
  rownames(sample_info_df) <- sample_info_df$index
  sample_info_df$index <- NULL # Remove redundant column
  meta_data <- sample_data(sample_info_df)

  # 3. Separate abundance data and transpose it
  otu_matrix <- raw_data_in %>%
    select(-index, -Group) %>%
    as.matrix()
  rownames(otu_matrix) <- rownames(sample_info_df)
  otu_tab <- otu_table(t(otu_matrix), taxa_are_rows = TRUE, errorIfNULL = TRUE)

  # 4. Parse taxonomy strings from column headers
  taxa_df <- data.frame(Taxon = colnames(otu_matrix)) %>%
    separate(Taxon,
             into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
             sep = ";",
             fill = "right") %>%
    mutate(across(everything(), ~ gsub("^[dkpcofgs]_[0-9]__", "", .))) %>%
    mutate(across(everything(), ~ na_if(., "")))
  rownames(taxa_df) <- colnames(otu_matrix)
  tax_tab <- phyloseq::tax_table(as.matrix(taxa_df))

  # Create an unrooted tree from taxa abundance profiles
  taxa_dist <- dist(otu_tab)
  phy_tree <- ape::nj(taxa_dist)

  # Build a phyloseq object
  physeq <- phyloseq(otu_tab, tax_tab, meta_data, phy_tree)

  # Filter taxa before running LEfSe
  physeq_ab_filtered <- prune_taxa(taxa_sums(physeq) > abundance_filter_threshold, physeq)
  physeq_filtered <- prune_taxa(rowSums(otu_table(physeq_ab_filtered) > 0) >= 3, physeq_ab_filtered)

  # Run LEfSe analysis
  lefse_result <- run_lefse(physeq_filtered,
    group = "Group",
    transform = transform,
    lda_cutoff = lda_cutoff,
    strict = strict,
    multigrp_strat = multigrp_strat,
    sample_min = sample_min_threshold,
    curv = curv
  )
  
  return(lefse_result)
}

### --- Execute Analysis --- ###
# Call the function to run the analysis and store the results
lefse_result <- run_lefse_analysis()

# Inspect the lefse results marker table:
marker <- print(marker_table(lefse_result))



