#'Generate Node Data Frame for Differential Cladogram
#'
#' This function performs differential abundance analysis using `MicrobiotaProcess::diff_analysis`
#' and prepares a data frame suitable for plotting a differential cladogram.
#' It identifies features (taxa) that are significantly different between specified groups.
#'
#' @param sub_physeq A phyloseq object, typically one that has been filtered and subsetted.
#' @param classgroup A character string specifying the name of the metadata column
#'                   in `sample_data(sub_physeq)` to use for grouping samples (e.g., "Group").
#' @param mlfun A character string specifying the machine learning function to use for
#'              differential analysis. Default is "lda" for LEfSe-like analysis.
#' @param ldascore A numeric value specifying the LDA score cutoff for significance.
#'                 Features with an absolute LDA score below this value will not be considered significant. Default is 3.
#' @param type A character string specifying the level of taxonomic resolution for the analysis
#'             (e.g., "species", "genus", "family"). Default is "species".
#' @param normalization A numeric value for normalizing abundance data. Default is 1,000,000.
#' @return A data frame containing the differential analysis results, including LDA scores,
#'         p-values, and enriched groups for each feature, ready for cladogram plotting.
crt_diffclad_obj <- function(obj, classgroup = "Group", mlfun = "lda", ldascore = 3,
                              type = "species", normalization = 1000000) {
  # 0. Check and load required packages
  # Ensure 'MicrobiotaProcess' is loaded for 'as.treedata'
  required_package <- c("MicrobiotaProcess", "dplyr")
  check_and_load_packages(required_package)

  # Perform differential abundance analysis using diff_analysis from MicrobiotaProcess
  # This function identifies features (taxa) that are significantly different
  # between the specified `classgroup` categories.
  message("Performing differential abundance analysis...")
  diffclad_obj <- diff_analysis(obj,
                                classgroup = classgroup,
                                mlfun = mlfun,
                                ldascore = ldascore,
                                normalization = normalization,
                                type = type,
                                clmin = 1
                                    )
  df_diffclad_obj <- as.data.frame(diffclad_obj) %>% rename(nodeLab = f) %>%
                      mutate(nodeLab = gsub("D_[0-9]__", "", nodeLab))
 

  # Return the result, which is a data frame containing the differential analysis
  # results, typically including LDA scores, p-values, and enriched groups for each feature.
  message("Done.")
  return(df_diffclad_obj) 
                               
}


