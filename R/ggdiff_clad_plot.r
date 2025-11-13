# Source required R scripts 
source("R/check_and_load_packages.r") # function for check and load packages
source("R/generate_physeq_obj.r") # function for generate phyloseq objects
source("R/diffclad_obj.r") # function for create diffclad input obj
source("R/diffclad_nodedf.r") # function for create diffclad input nodedf


# Load required packages (The sourced R scripts define its own required packages)
# The following required packages is only for this R script.
required_pkgs <- c("MicrobiotaProcess", "phyloseq", "ggplot2", "ggtree")
check_and_load_packages(required_pkgs) # Load necessary packages, replacing ggtext with gridExtra


# Generate the initial phyloseq object from raw data
physeq_obj <- generate_physeq_obj(raw_data = "data/processed/ALL-update.csv",
                                  META_cols = c("index", "Group"),
                                  sample_id_col = "index")

# Subset the phyloseq object based on Group, abundance and prevalence thresholds.
# The subsetting parameters are defined in generate_physeq_obj.r

message("Filtering physeq's taxa on abundance_filter_threshold...")
# First, filter taxa based on abundance across all samples
physeq_ab_filtered <- prune_taxa(taxa_sums(physeq_obj) > 600, physeq_obj)
  
message("Filtering physeq's taxa on prevalence_filter_threshold...")
# Second, filter out taxa that are not present in at least 3 samples across the full dataset.
physeq_prev_filtered <- prune_taxa(rowSums(otu_table(physeq_ab_filtered) > 0) >= 3, physeq_ab_filtered)

message("Subsetting physeq...")
# Finally, subset the fully filtered data to the groups of interest for the pairwise comparison
# We use .data[[class_col]] to correctly evaluate the column name passed as a string
  
subset_physeq <- phyloseq::subset_samples(physeq_prev_filtered, Group %in% c("EHI", "EHI.abt"))
message("Done")

# Create the ggdiffclade input (nodedf)
diffclad_nodedf <- crt_diffclad_nodedf(obj = subset_physeq,classgroup = "Group", mlfun = "lda", ldascore = 3,
                              type = "species", normalization = 1000000)

# Create the ggdiffclade input (obj)
diffclad_obj <- crt_diffclad_obj(subset_physeq = subset_physeq)


# Create the ggdiffclade plot

p <- ggdiffclade(obj = diffclad_obj,
                nodedf = diffclad_nodedf,
                factorName = "Group",
                layout = "circular",
                skpointsize = 0.6,
                cladetext = 2,
                linewd = 0.2,
                reduce = TRUE,
                taxlevel = 1) +
    scale_fill_manual(values = c("#00AED7", "#009E73")) +
    guides(color = guide_legend(
        keywidth = 0.1,
        keyheight = 0.6,
        order = 3,
        ncol = 1
    )) +
    theme(
        panel.background = element_rect(fill = NA),
        legend.position = c(1.15, 0.55),
        plot.margin = margin(0, 4.5, 0, 0, "cm"),
        legend.spacing.y = unit(0.5, "cm"),
        legend.title = element_text(size = 9.5),
        legend.text = element_text(size = 7.5),
        legend.box.spacing = unit(0, "cm"),
        plot.caption = element_text(
            size = 14,
            face = "bold",
            hjust = .7)
        ) +     
    labs(caption = "Differential Aboundance Cladogram - EHI vs EHI.abt")


print(p)

ggsave("reports/figures/diff_clade_EHI.png", plot = p, width = 12, height = 12, dpi = 600)
