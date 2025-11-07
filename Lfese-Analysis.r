### Load Required Libraries ###

# This script assumes all required packages have been installed.
# Use install.packages("package_name") or BiocManager::install("package_name") if any are missing.

library(phyloseq)
library(microbiomeMarker) # For LEfSe analysis
library(dplyr)
library(data.table)       # For the fread function
library(ape)

### Create the phyloseq object from the CSV file
# This block reads the raw data, parses it into abundance, taxonomy, and metadata,
# and constructs the phyloseq object required for all subsequent analysis.

# 1. Read the data using data.table for speed
raw_data_in <- fread("ALL-update.csv")

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
taxa_names_vec <- colnames(otu_matrix)
tax_matrix_list <- strsplit(taxa_names_vec, ";")

# Clean prefixes (e.g., "D_0__", "k__") and create a matrix
tax_matrix_cleaned <- sapply(tax_matrix_list, function(x) {
  gsub("^[dkpcofgs]__", "", x)
})
tax_matrix_cleaned <- t(tax_matrix_cleaned)
tax_matrix_cleaned[tax_matrix_cleaned == ""] <- NA

# Assign standard rank names
rank_names <- c(
  "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"
)
colnames(tax_matrix_cleaned) <- rank_names[seq_len(ncol(tax_matrix_cleaned))]
rownames(tax_matrix_cleaned) <- taxa_names_vec

# Create the tax_table for the phyloseq object
tax_tab <- phyloseq::tax_table(tax_matrix_cleaned)

# Create an unrooted tree from taxa abundance profiles
taxa_dist <- dist(otu_tab)
phy_tree <- ape::nj(taxa_dist)

### Build a phyloseq object
physeq <- phyloseq(otu_tab, tax_tab, meta_data, phy_tree)

### Run LEfSe analysis by using run_lefse
lefse_result <- run_lefse(physeq, group = "Group", lda_cutoff = 3, strict = "1")


### Plot
library(ggtree)
library(ggtreeExtra) # To support the outer rings of a tree in a circular layout
library(ggplot2)

# Prepare data for plotting
Psq <- physeq
Psq <- prune_taxa(taxa_sums(Psq) > 1000, Psq)

# Fix: merge_samples() is fundamentally broken when sample_data contains non-numeric columns
# Instead of using merge_samples, we manually aggregate OTU table and create new phyloseq object
# This avoids the "NAs introduced by coercion" warnings and dimension corruption issues

# Get the grouping variable
grouping_var <- as.character(sample_data(Psq)$Group)

# Manually aggregate OTU table by group
otu_mat <- as(otu_table(Psq), "matrix")
if (!taxa_are_rows(Psq)) {
  otu_mat <- t(otu_mat)
}

# Sum abundances within each group
groups <- unique(grouping_var)
merged_otu <- matrix(0, nrow = nrow(otu_mat), ncol = length(groups))
rownames(merged_otu) <- rownames(otu_mat)
colnames(merged_otu) <- groups

for (grp in groups) {
  samples_in_group <- which(grouping_var == grp)
  if (length(samples_in_group) == 1) {
    merged_otu[, grp] <- otu_mat[, samples_in_group]
  } else {
    merged_otu[, grp] <- rowSums(otu_mat[, samples_in_group])
  }
}

# Create new phyloseq object with merged data
merged_otu_table <- otu_table(merged_otu, taxa_are_rows = TRUE)
merged_sample_data <- data.frame(
  Group = groups,
  row.names = groups,
  stringsAsFactors = FALSE
)
merged_sample_data <- sample_data(merged_sample_data)

# Build the merged phyloseq object (reuse tax_table and tree from original)
MergedPsq <- phyloseq(
  merged_otu_table,
  tax_table(Psq),
  merged_sample_data,
  phy_tree(Psq)
)

# Rarefy merged samples to normalize group-level abundances for visualization
MergedPsq <- rarefy_even_depth(MergedPsq, rngseed = 394582)

# Create melted data frame with corrected column names
melt_simple <- psmelt(MergedPsq) %>% 
  filter(Abundance < 120) %>% 
  select(OTU, val = Abundance)

# Create fan tree box plot with phyla coloring
# Extract the tree from phyloseq object before merging (to preserve tree structure)
tree <- phy_tree(Psq)

# Get taxonomy information for coloring by Phylum
# Create a data frame with taxonomy info where label column matches tree tip labels
tax_df <- as.data.frame(tax_table(Psq))
# Add label column that exactly matches tree$tip.label for proper joining
tax_df$label <- rownames(tax_df)

# Debug: Print some diagnostic information
cat("Number of tree tips:", length(tree$tip.label), "\n")
cat("Number of taxa:", nrow(tax_df), "\n")
cat("First 5 tree tip labels:", head(tree$tip.label, 5), "\n")
cat("First 5 taxonomy labels:", head(tax_df$label, 5), "\n")
cat("Unique Phyla found:", length(unique(tax_df$Phylum[!is.na(tax_df$Phylum)])), "\n")
cat("Phylum values sample:", head(unique(tax_df$Phylum), 10), "\n")

# Create the base ggtree plot with fan layout
# Using layout = "fan" with open.angle = 10 for a wedge-shaped tree
p <- ggtree(tree, layout = "fan", open.angle = 10)

# Join the taxonomy data - ggtree uses the 'label' column to match tree tip labels
p <- p %<+% tax_df

# Add colored tip points by Phylum
p <- p + geom_tippoint(aes(color = Phylum), size = 2, alpha = 0.8) +
  theme(legend.position = "right",
        legend.title = element_text(size = 11, face = "bold"),
        legend.text = element_text(size = 9),
        legend.key.size = unit(0.5, "cm")) +
  guides(color = guide_legend(override.aes = list(size = 4, alpha = 1)))

# Add abundance data as boxplot layer using ggtreeExtra
p <- p + 
  geom_fruit(
    data = melt_simple,
    geom = geom_boxplot,
    mapping = aes(y = OTU, x = val),
    orientation = "y",
    pwidth = 0.3,
    axis.params = list(
      axis = "x",
      text.size = 3,
      title = "Abundance"
    )
  )

# Display the plot
print(p)

# Optionally save the plot
# ggsave("fan_tree_boxplot.pdf", plot = p, width = 12, height = 12)