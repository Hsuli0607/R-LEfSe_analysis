#' Generate Phyloseq Object
#'
#' This function encapsulates the process from reading input data
#' (a CSV file) to generating a phyloseq object.
#'
#' @param raw_data Path to the input CSV file. The file should have
#'                 samples as rows.
#' @param META_cols A character vector specifying the names of all metadata
#'                  columns (e.g., "index", "Group", "Age").
#' @param sample_id_col A single character string specifying which of the
#'                      `META_cols` contains the unique sample IDs. This
#'                      column will be used to set the rownames.
#'
#' @return A phyloseq object containing the OTU table, sample data, and
#'         taxonomy table. No phylogenetic tree is included.
#' @export
#'
#' @examples
#' \dontrun{
#' # Assuming 'ALL-update.csv' is in a 'data/processed' folder
#' # and has columns "index", "Group", and then many taxa columns.
#' pseq <- generate_physeq_obj(
#'   raw_data = "data/processed/ALL-update.csv",
#'   META_cols = c("index", "Group"),
#'   sample_id_col = "index"
#' )
#' }
generate_physeq_obj <- function(
    raw_data = "data/processed/ALL-update.csv",
    META_cols = c("index", "Group"),
    sample_id_col = "index") {

  # 0. Check and load required packages
  # Ensure you have a function `check_and_load_packages` or load manually
  required_pkgs <- c("data.table", "dplyr", "tidyr", "phyloseq", "magrittr")
  check_and_load_packages(required_pkgs)

  # Ensure sample_id_col is one of the META_cols
  if (!sample_id_col %in% META_cols) {
    stop("`sample_id_col` must be one of the names in `META_cols`.")
  }

  # Example of manual loading if check_and_load_packages doesn't exist:
  # invisible(lapply(required_pkgs, function(pkg) {
  #   if (!require(pkg, character.only = TRUE)) {
  #     warning(paste("Package", pkg, "not found. Please install it."))
  #   }
  # }))
  
  message("Loading data...")
  # 1. Read the data using data.table for speed
  raw_data_in <- data.table::fread(raw_data)
  
  message("Processing sample metadata...")
  # 2. Separate sample metadata
  sample_info_df <- raw_data_in %>%
    dplyr::select(dplyr::all_of(META_cols)) %>%
    as.data.frame()
  
  # Set rownames from the specified sample ID column
  rownames(sample_info_df) <- sample_info_df[[sample_id_col]]
  
  # Remove the redundant sample ID column (as pointed out)
  sample_info_df[[sample_id_col]] <- NULL
  
  # Create the phyloseq sample_data object
  meta_data <- phyloseq::sample_data(sample_info_df)
  
  message("Processing abundance data (OTU table)...")
  # 3. Separate abundance data
  otu_matrix <- raw_data_in %>%
    dplyr::select(-dplyr::all_of(META_cols)) %>%
    as.matrix()
  
  # Set sample names (rows)
  rownames(otu_matrix) <- rownames(sample_info_df)
  
  # Transpose and create the otu_table
  # Taxa should be rows in the final otu_table
  otu_tab <- phyloseq::otu_table(t(otu_matrix), taxa_are_rows = TRUE)
  
  message("Processing taxonomy...")
  # 4. Parse taxonomy strings from column headers
  taxa_df <- data.frame(Taxon = colnames(otu_matrix)) %>%
    tidyr::separate(
      Taxon,
      into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
      sep = ";",
      fill = "right", # Fill missing ranks with NA
      remove = TRUE    # Remove the original 'Taxon' column
    ) %>%
    # Clean up prefixes like k__, p__, c__ etc.
    dplyr::mutate(dplyr::across(everything(), ~ gsub("^[dkpcofgs]__", "", .))) %>%
    # Convert empty strings "" to NA
    dplyr::mutate(dplyr::across(everything(), ~ dplyr::na_if(., "")))
  
  # Set taxon names (rows)
  rownames(taxa_df) <- colnames(otu_matrix)
  
  # Create the phyloseq tax_table
  tax_tab <- phyloseq::tax_table(as.matrix(taxa_df))
  
  message("Building phyloseq object...")
  # 5. Build a phyloseq object
  # Removed the phy_tree generation as it was not a valid phylogenetic tree.
  physeq <- phyloseq::phyloseq(otu_tab, tax_tab, meta_data)
  
  message("Done.")
  return(physeq)
}

 