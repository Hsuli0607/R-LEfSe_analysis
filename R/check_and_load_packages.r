#' Check and Load Required Packages Function (Original Version)
#'
#' This function checks if specified R packages are installed. If any are not
#' installed, it stops execution and prompts the user to install them.
#' Otherwise, it loads all specified packages, suppressing startup messages.
#'
#' @param pkgs A character vector of package names to check and load.
#' @return Invisible NULL. Stops execution if packages are missing.
check_and_load_packages <- function(pkgs) {
  installed_pkgs <- rownames(installed.packages())
  missing_pkgs <- pkgs[!pkgs %in% installed_pkgs]

  if (length(missing_pkgs) > 0) {
    message("Error: The following required packages are not installed:")
    message(paste(missing_pkgs, collapse = ", "))
    stop("Please install the missing packages before proceeding.", call. = FALSE)
  }

  suppressPackageStartupMessages({
    lapply(pkgs, library, character.only = TRUE)
  })
  invisible(NULL)
}