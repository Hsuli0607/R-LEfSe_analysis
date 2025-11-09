### Load Required Libraries ###
# This script assumes all required packages are installed.
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

# This function takes a path to a CSV file containing microbiome data,
# performs LEfSe analysis, and generates a fan tree plot with barplots.

# @param input_path Path to the input CSV file. 
# The CSV should have an 'index' column with sample names, a 'Group' column for sample groups,
# and the rest of the columns as taxa with their abundances.
# @param lda_cutoff LDA score cutoff for LEfSe analysis. Default is 3.
# @param strict Strictness of the LEfSe analysis. Default is "1".
# @param rarefy_seed Seed for reproducibility of rarefaction. Default is 394582.
# @param transform The transformation to apply to the data before analysis. Default is "log10p".
# @param abundance_filter_threshold Threshold for filtering taxa by abundance. Default is 1000.
# return A ggtree plot object P.

# Instructioins how to run this function:
# Simply call it with the desired parameters. 
# plot <- run_lefse_analysis_and_plot(input_path = "data/processed/ALL-update.csv")
# print(plot)
 
run_lefse_analysis_and_plot <- function(input_path = "ALL-update.csv",
                                        lda_cutoff = 3,
                                        strict = "1",
                                        rarefy_seed = 394582,
                                        abundance_filter_threshold = 1000,
                                        transform = "log10p",
                                        multigrp_strat = FALSE,
                                        curv = FALSE,
                                        sample_min_threshold = 3) {

  ### Create the phyloseq object from the CSV file
  # This block reads the raw data, parses it into abundance, taxonomy, and metadata,
  # and constructs the phyloseq object required for subsequent analysis.
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

  # Explicitly remove branch lengths as they are not used in this visualization
  phy_tree$edge.length <- NULL

  # Explicitly remove branch lengths to prevent ggtree from drawing its own scale
  phy_tree$edge.length <- NULL
  ### Build a phyloseq object
  physeq <- phyloseq(otu_tab, tax_tab, meta_data, phy_tree)

  # Filter taxa before running LEfSe to avoid issues with low-abundance features.
  physeq_filtered <- prune_taxa(taxa_sums(physeq) > abundance_filter_threshold, physeq)

  ### Run LEfSe analysis by using run_lefse
  lefse_result <- run_lefse(physeq_filtered,
    group = "Group",
    transform = transform,
    lda_cutoff = lda_cutoff,
    strict = strict,
    multigrp_strat = FALSE
  )
  
  # Extract the marker data to link taxa to their enriched group
  marker_data <- data.frame(marker_table(lefse_result))
  marker_data_enrich <- marker_data %>%
    select(OTU = feature, enrich_group)

  ### Plot
  # Prepare data for plotting the tree and associated bar plots
  Psq <- physeq_filtered
  tree <- phy_tree(Psq)
  taxdf <- as.data.frame(tax_table(Psq))
  taxdf$OTU <- rownames(taxdf) # Use rownames as the OTU identifier


  # Create a single, comprehensive data frame for plotting
  plot_data <- psmelt(Psq) %>% # Using psmelt to get a long-form data frame
    group_by(OTU) %>%
    summarise(mean_abundance = mean(log10(Abundance + 1)), .groups = 'drop') %>%
    left_join(taxdf, by = "OTU") %>%
    # Clean up Phylum names for the legend
    mutate(Phylum = if_else(is.na(Phylum) | Phylum %in% c("_", "D_1__", "__"), "Unknown", Phylum)) %>%
    # Add the enrichment group information for significant taxa
    left_join(marker_data_enrich, by = "OTU")

  # Build ggtree with a fan layout and add taxon data
  p <- ggtree(tree, layout = "fan", open.angle = 20, branch.length = "none") %<+% plot_data +
    geom_tippoint(aes(color = Phylum), size = 1.5, show.legend = FALSE)
  # Start with a completely blank theme and add only what is needed
  p <- p + theme_void() +
    theme(
      panel.background = element_rect(fill = "beige", colour = "beige"),
      legend.position = "right",
      legend.title = element_text(size = 11, face = "bold"),
      legend.text = element_text(size = 9),
      legend.key.size = unit(0.5, "cm")
    )

  # Add the abundance bar plot, colored by Phylum
  p <- p +
    geom_fruit(
      # No 'data' argument needed; it inherits from the main plot
      geom = geom_col,
      mapping = aes(
        y = OTU,
        x = mean_abundance,
        fill = Phylum
      ),
      width = 0.38,
      orientation = "y",
      axis.params = list(
        axis = "x",
        text.size = 2,
        title = "Mean Abundance (log10)",
        vjust = 1
      ),
      grid.params = list(linetype = "dotted", color = "grey80") # Add grid lines
    )

  # Add a second layer to show the enrichment group for each significant taxon
  p <- p +
    geom_fruit(
      geom = geom_tile,
      mapping = aes(
        y = OTU,
        fill = enrich_group # Color tiles by the enrichment group
      ),
      width = 0.1,
      offset = 0.08 # Adjust offset to place it next to the bar plot
    ) +
    # Add a new color scale for the enrichment groups
    ggnewscale::new_scale_fill() +
    scale_fill_viridis_d(name = "Enriched Group", na.value = "white", option = "C")

  # Add final theming and a single, styled legend for Phylum
  p <- p +
    labs(fill = "Phyla") + # Set the legend title
    guides(fill = guide_legend(ncol = 1)) # Ensure single-column legend
    return(p)
}
### Run the function and print and save plot:
plot <- run_lefse_analysis_and_plot(input_path = "data/processed/ALL-update.csv")
print(plot)
ggsave("reports/figures/plot.png", plot = plot, width = 12, height = 12, dpi = 300)git branch -vv