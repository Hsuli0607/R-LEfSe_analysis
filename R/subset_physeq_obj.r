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

subset_physeq_obj <- function(
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