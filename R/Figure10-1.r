# Fixed script to plot a fan tree with boxplot "fruits" using GlobalPatterns
# - Ensures tip metadata has an explicit 'label' column before using %<+%
# - Ensures geom_fruit data uses 'label' to match tree tip labels
# - Minor safety: sets seed for rarefaction and avoids factor surprises

library(ggtreeExtra)
library(ggtree)
library(phyloseq)
library(dplyr)
library(ggplot2)

# Load data
data("GlobalPatterns")
GP <- GlobalPatterns

# Filter taxa, add a helper sample variable
GP <- prune_taxa(taxa_sums(GP) > 600, GP)
sample_data(GP)$human <- get_variable(GP, "SampleType") %in% c("Feces", "Skin")

# Merge samples by SampleType, rarefy, and agglomerate at Order level
mergedGP <- merge_samples(GP, "SampleType")
set.seed(394582)
mergedGP <- rarefy_even_depth(mergedGP, rngseed = 394582, verbose = FALSE)
mergedGP <- tax_glom(mergedGP, taxrank = "Order")

# Extract tree
tree <- phy_tree(mergedGP)

# Build a clean taxon metadata data.frame and ensure there's a 'label' column
taxdf <- as.data.frame(tax_table(mergedGP), stringsAsFactors = FALSE)
taxdf$label <- taxa_names(mergedGP)   # explicit tip label column
# Optionally reorder so 'label' is first (not required, but clearer)
taxdf <- taxdf[, c("label", setdiff(names(taxdf), "label"))]

# Prepare the psmelt data used for the geom_fruit boxplots
melt_df <- psmelt(mergedGP)

# Inspect these if you want to confirm:
# head(taxa_names(mergedGP))
# head(unique(melt_df$OTU))

# In psmelt after merge_samples, OTU contains the taxa names that should match taxa_names(mergedGP).
# Rename OTU to label so geom_fruit can match to tree tip labels.
# Note: Don't include Phylum and Order here since they're already in tree data via %<+% taxdf
melt_simple <- melt_df %>%
  filter(Abundance < 120) %>%
  select(label = OTU, Sample = Sample, val = Abundance)

# If your grouping for boxplots is by SampleType, use Sample (which is the merged SampleType).
# If you prefer grouping by something else, change group mapping below.

# Build ggtree and attach taxon metadata with %<+%
p <- ggtree(tree, layout = "fan", open.angle = 10) %<+% taxdf +
  geom_tippoint(aes(color = Phylum), size = 1.5)

# Add the boxplot "fruit". Use data that contains 'label' column matching tip labels.
p <- p +
  geom_fruit(
    data = melt_simple,
    geom = geom_boxplot,
    mapping = aes(
      y = label,
      x = val,
      group = label, # Group by taxon to get one boxplot per tip
      fill = Phylum  # Use Phylum for fill color
    ),
    linewidth = 0.2,
    outlier.size = 0.5,
    outlier.stroke = 0.08,
    outlier.shape = 21,
    axis.params = list(
      axis = "x",
      text.size = 1.8,
      hjust = 1,
      vjust = 0.5,
      nbreak = 3
    ),
    grid.params = list()
  )

# Final theming
p <- p +
  # Use the same aesthetic for color and fill to unify the legend
  labs(color = "Phyla", fill = "Phyla") +
  guides(fill = guide_legend(ncol = 1), color = guide_legend(ncol = 1)) +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 11, face = "bold"),
    legend.text = element_text(size = 9),
    legend.key.size = unit(0.5, "cm"),
    panel.background = element_rect(fill = "beige", colour = "beige"),
    panel.grid.major = element_line(colour = "grey90")
  )

# Print plot
print(p)

# save plot to 
ggsave("reports/figures/Figure-10.png",
  plot = p, width = 12, height = 12, dpi = 300
)
 
