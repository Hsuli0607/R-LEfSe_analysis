### Load Required Libraries ###
# This script assumes all required packages have been installed.
# Use install.packages("package_name") or BiocManager::install("package_name") if any are missing.

library(phyloseq)
library(microbiomeMarker) # For LEfSe analysis
library(dplyr)
library(data.table)       # For the fread function
library(tidyr)
library(ape)



### Generate Phyloseq Object Function

#' This function encapsulates the from reading 
#' the input data to generate phyloseq object.
#' @param raw_data Path to the input CSV file.
#' @Return return Phyloseq object

run_generate_physeq <- function(
                                 raw_data= "data/processed/ALL-update.csv"
                                 ) {
  required_pkgs <- c("data.table", "dplyr", "tidyr", "phyloseq", "ape")
  check_and_load_packages(required_pkgs)
  # 1. Read the data using data.table for speed
  raw_data_in <- fread(raw_data)

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

  return(physeq)
}

physeq <- run_generate_physeq()

### Subset And Filtering Phyloseq_Object Function
#'
#' This function first filters taxa by total abundance and prevalence
#' across the *entire* dataset, and *then* subsets the samples
#' based on a specific metadata column and values.
#'
#' @param physeq_object The input phyloseq object.
#' @param class_col A string (character) specifying the name of the 
#'                  metadata column to use for subsetting (e.g., "Group").
#' @param subset_vals A vector of values to keep from the `class_col` 
#'                    (e.g., c("EHI", "EHI.abt")).
#' @param abundance_filter_threshold The minimum total abundance a taxon 
#'                                 must have across all samples to be kept.
#' @param prevalence_filter_threshold The minimum number of samples a taxon
#'                                    must be present in (count > 0) to be kept.
#'
#' @return A new, filtered, and subsetted phyloseq object.

physeq_subset <- function(
                           physeq_object = run_generate_physeq(),
                           class_col = "Group",
                           subset_vals = c("EHI", "EHI.abt"),
                           abundance_filter_threshold = 1000,
                           prevalence_filter_threshold = 3
                         ) {
  required_pkgs <- c("phyloseq")
  check_and_load_packages(required_pkgs)
  
  # First, filter taxa based on abundance across all samples
  physeq_ab_filtered <- prune_taxa(taxa_sums(physeq) > abundance_filter_threshold, physeq)
  
  # Second, filter out taxa that are not present in at least 3 samples across the full dataset.
  physeq_prev_filtered <- prune_taxa(rowSums(otu_table(physeq_ab_filtered) > 0) >= prevalence_filter_threshold, physeq_ab_filtered)

  # Finally, subset the fully filtered data to the groups of interest for the pairwise comparison
  physeq_filtered <- subset_samples(physeq_prev_filtered, class_col %in% subset_vals)
  
  return(physeq_filtered)
}

physeq_filtered <- physeq_subset(physeq)

#' LEfSe Analysis Function
#'
#' @param lda_cutoff LDA score cutoff for LEfSe analysis. Default is 3.
#' @param strict Strictness of the LEfSe analysis. Default is "1".
#' @param transform The transformation to apply to the data. Default is "log10p".
#' @param sample_min_threshold Minimum number of samples a taxon must be present in. Default is 3.
#' @param curv Whether to use a curved LDA model. Default is TRUE.
#' @return A `microbiomeMarker` object containing the LEfSe results.

### LEfSe Analysis Function
lefse_analysis <- function(
                            physeq_object = physeq_filtered,
                            transform = "log10p",
                            lda_cutoff = 3,
                            strict = "1",
                            multigrp_strat = FALSE,
                            sample_min_threshold = 3,
                            curv = FALSE
                        ) {
  required_pkgs <- c("microbiomeMarker")
  check_and_load_packages(required_pkgs)
  lefse_result <- run_lefse(physeq_object,
                             transform = transform,
                             lda_cutoff = lda_cutoff,
                             strict = strict,
                             multigrp_strat = multigrp_strat,
                             sample_min = sample_min_threshold,
                             curv = curv
                           )
return(lefse_result)
}

lefse_result <- lefse_analysis(physeq_object = physeq)





### Plot Diff_Cladogram Function
