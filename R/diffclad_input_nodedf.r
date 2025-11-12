#' Generate Diff_Cladogram Input Object: obj
#'#' This function generates a `treedata` object from a phyloseq object's
#' taxonomy table, which can then be used as input for `ggdiffclade`.
#'
#' @param subset_physeq A phyloseq object, typically one that has been
#'                      filtered and subsetted.
#'
#' @return A `treedata` object derived from the taxonomy table of the
#'         input phyloseq object.
#' @export
#'
#' @examples
#' \dontrun{
#' # Assuming 'physeq_filtered' is a phyloseq object
#' # diffclade_obj <- ggdiffclade_input_obj(subset_physeq = physeq_filtered)
#' }
crt_diffclad_nodedf <- function(subset_physeq)
                                  {
  # 0. Check and load required packages
  # Ensure 'MicrobiotaProcess' is loaded for 'as.treedata'
  required_package <- c("MicrobiotaProcess", "phyloseq", "stringr")
  check_and_load_packages(required_package)

  # Extract the taxonomy table from the phyloseq object
  message("Extracting taxonomy table...")
  tax_tab <- phyloseq::tax_table(subset_physeq)

  # Convert the taxonomy table to dataframe
  message("Converting to dataframe...")
  df_tax_tab <- as.data.frame(tax_tab)

  # Transform to nodedf
  message("Converting to nodedf...")

  tax_cols <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
  diffclad_nodedf <- df_tax_tab %>%
  mutate(across(all_of(tax_cols), ~ str_remove(., "^.*__")))

  message("Done.")
  return(diffclad_nodedf)
 
}

