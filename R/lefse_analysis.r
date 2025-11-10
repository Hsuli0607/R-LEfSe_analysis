#' LEfSe Analysis Function
#'
#' @param lda_cutoff LDA score cutoff for LEfSe analysis. Default is 3.
#' @param strict Strictness of the LEfSe analysis. Default is "1".
#' @param transform The transformation to apply to the data. Default is "log10p".
#' @param sample_min_threshold Minimum number of samples a taxon must be present in. Default is 3.
#' @param curv Whether to use a curved LDA model. Default is TRUE.
#' @return A `microbiomeMarker` object containing the LEfSe results.

### LEfSe Analysis Function
lefse_analysis <- function(
                            physeq_object = physeq_filtered,
                            transform = "log10p",
                            lda_cutoff = 3,
                            strict = "1",
                            multigrp_strat = FALSE,
                            sample_min_threshold = 3,
                            curv = FALSE
                        ) {
  required_pkgs <- c("microbiomeMarker")
  check_and_load_packages(required_pkgs)
  lefse_result <- run_lefse(physeq_object,
                             transform = transform,
                             lda_cutoff = lda_cutoff,
                             strict = strict,
                             multigrp_strat = multigrp_strat,
                             sample_min = sample_min_threshold,
                             curv = curv
                           )
return(lefse_result)
}