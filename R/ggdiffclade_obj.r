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


ggdiffclade_obj <- function(subset_physeq = subset_physeq)
                                  {
  # 0. Check and load required packages
  # Ensure 'MicrobiotaProcess' is loaded for 'as.treedata'
  required_package <- c("MicrobiotaProcess")
  check_and_load_packages(required_package)

  # Extract the taxonomy table from the phyloseq object

Tax_tab <- tax_table(subset_physeq)
ggdiffclad_obj <- as.treedata(Tax_tab)

return(ggdiffclade_obj)

}