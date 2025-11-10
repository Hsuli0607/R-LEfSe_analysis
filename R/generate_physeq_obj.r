## Generate Phyloseq Object Function

#' This function encapsulates the from reading 
#' the input data to generate phyloseq object.
#' @param raw_data Path to the input CSV file.
#' @Return return Phyloseq object

generate_physeq_obj <- function(
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