# File updated to resolve potential environment caching issues.
### Load Required Libraries ###
# This script assumes all required packages have been installed.
# Use install.packages("package_name") or BiocManager::install("package_name") if any are missing.
library(phyloseq)
library(microbiomeMarker) # For LEfSe analysis
library(dplyr)
library(tidyr)
library(data.table)       # For the fread function
library(ape)
library(ggtree)
library(ggtreeExtra) # To support the outer rings of a tree in a circular layout
library(ggplot2)
#' Run LEfSe Analysis and Generate a Plot
#'
#' This function takes a path to a CSV file containing microbiome data,
#' performs LEfSe analysis, and generates a fan tree plot with boxplots.
#'
#' @param input_path Path to the input CSV file. The CSV should have an 'index' column
#'   with sample names, a 'Group' column for sample groups, and the rest of the
#'   columns as taxa with their abundances.
#' @param lda_cutoff LDA score cutoff for LEfSe analysis. Default is 3.
#' @param strict Strictness of the LEfSe analysis. Default is "1".
#' @param rarefy_seed Seed for reproducibility of rarefaction. Default is 394582.
#' @param abundance_filter_threshold Threshold for filtering taxa by abundance. Default is 1000.
#'
#' @return A ggtree plot object.
#' @export
run_lefse_analysis_and_plot <- function(input_path = "ALL-update.csv",
                                        lda_cutoff = 3,
                                        strict = "1",
                                        rarefy_seed = 394582,
                                        abundance_filter_threshold = 1000) {

  ### Create the phyloseq object from the CSV file
  # This block reads the raw data, parses it into abundance, taxonomy, and metadata,
  # and constructs the phyloseq object required for all subsequent analysis.
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
  # Phyloseq requires taxa as rows and samples as columns
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
    mutate(across(everything(), ~ gsub("^[dkpcofgs]__", "", .))) %>%
    mutate(across(everything(), ~ na_if(., "")))
  rownames(taxa_df) <- colnames(otu_matrix)
  tax_tab <- phyloseq::tax_table(as.matrix(taxa_df))

  # Create an unrooted tree from taxa abundance profiles
  taxa_dist <- dist(otu_tab)
  phy_tree <- ape::nj(taxa_dist)
  ### Build a phyloseq object
  physeq <- phyloseq(otu_tab, tax_tab, meta_data, phy_tree)

  # Filter taxa before running LEfSe to avoid issues with low-abundance features.
  physeq_filtered <- prune_taxa(taxa_sums(physeq) > abundance_filter_threshold, physeq)

  ### Run LEfSe analysis by using run_lefse
  lefse_result <- run_lefse(physeq_filtered,
    group = "Group", lda_cutoff = lda_cutoff, strict = strict,
    multigrp_strat = TRUE
  )

  ### Plot
  # Prepare data for plotting
  Psq <- physeq_filtered
  tree <- phy_tree(Psq)
  taxdf <- as.data.frame(tax_table(Psq))

  # Calculate mean abundance for the bar plot layer
  abundance_data <- psmelt(Psq) %>%
    group_by(OTU) %>%
    summarise(mean_abundance = mean(Abundance))

  # build ggtree
  p <- ggtree(tree, layout = "fan", open.angle = 10) %<+% taxdf +
    geom_tippoint(aes(color = Phylum), size = 1.5) +
    theme(
      legend.position = "right",
      legend.title = element_text(size = 9, face = "bold"),
      legend.text = element_text(size = 9),
      legend.key.size = unit(0.5, "cm"),
      panel.background = element_rect(fill = "#F5F5DC", colour = "#F5F5DC"), # Beige color
      panel.grid.major = element_line(colour = "grey90")
    )

  # add the abundance boxplot
  p <- p +
    geom_fruit(
      data = abundance_data,
      geom = geom_col,
      mapping = aes(
        y = OTU,
        x = mean_abundance
      ),
      fill = "grey20", # Set a static fill color for the bars
      orientation = "y",
      pwidth = 0.3,
      axis.params = list(
        axis = "x",
        text.size = 2,
        title = "Mean Abundance"
      ),
      grid.params = list()
    )

  # Add final theming and legend control
  p <- p +
    labs(color = "Phylum") +
    guides(color = guide_legend(override.aes = list(size = 3), ncol = 1))

  return(p)
}
# Example of how to run the function:
plot_result <- run_lefse_analysis_and_plot(input_path = "data/processed/ALL-update.csv")
print(plot_result)
ggsave("reports/figures/plot.png", plot = plot_result, width = 12, height = 12, dpi = 300)