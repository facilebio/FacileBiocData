#' @include api.R
NULL

#' @export
#' @importFrom utils packageVersion
#' @importClassesFrom MultiAssayExperiment MultiAssayExperiment
setClass("FacileMultiAssayExperiment",
         contains = c("FacileBiocDataStore", "MultiAssayExperiment"))

#' @export
#' @noRd
facilitate.MultiAssayExperiment <- function(x, ...) {
  reqpkg("MultiAssayExperiment")
}

# bioc data retrieval methods --------------------------------------------------

#' @noRd
fdata.MultiAssayExperiment <- function(x, assay_name = default_assay(x), ...) {
  reqpkg("MultiAssayExperiment")
  stop("MultiAssayExpermient support not yet implemented")
}

#' @noRd
pdata.MultiAssayExperiment <- function(x, ...) {
  reqpkg("MultiAssayExperiment")
  stop("MultiAssayExpermient support not yet implemented")
}

#' @noRd
adata.MultiAssayExperiment <- function(x, name = default_assay(x), ...) {
  reqpkg("MultiAssayExperiment")
  stop("MultiAssayExpermient support not yet implemented")
}

# facile -----------------------------------------------------------------------

#' @noRd
#' @export
assay_names.MultiAssayExperiment <- function(x, ...) {
  reqpkg("MultiAssayExperiment")
  stop("MultiAssayExpermient support not yet implemented")
}
