#' @include api.R
NULL

#' @export
#' @importFrom utils packageVersion
#' @importClassesFrom MultiAssayExperiment MultiAssayExperiment
#' @importClassesFrom S4Vectors DataFrame SimpleList
setClass("FacileMultiAssayExperiment",
         slots = c(facile = "list"),
         contains = c("FacileBiocDataStore", "MultiAssayExperiment"),
         prototype = prototype(facile = list()))

#' @export
#' @noRd
#' @rdname facilitate
#' @method facilitate MultiAssayExperiment
facilitate.MultiAssayExperiment <- function(x, ...) {

}

# bioc data retrieval methods --------------------------------------------------

#' @noRd
fdata.MultiAssayExperiment <- function(x, ...) {
  if (!requireNamespace("MultiAssayExperiment", quietly = TRUE)) {
    stop("MultiAssayExperiment package required, please install it.",
         call. = FALSE)
  }
  stop("MultiAssayExpermient support not yet implemented")
}

#' @noRd
pdata.MultiAssayExperiment <- function(x, ...) {
  if (!requireNamespace("MultiAssayExperiment", quietly = TRUE)) {
    stop("MultiAssayExperiment package required, please install it.",
         call. = FALSE)
  }
  stop("MultiAssayExpermient support not yet implemented")
}

#' @noRd
adata.MultiAssayExperiment <- function(x, name = NULL, ...) {
  if (!requireNamespace("MultiAssayExperiment", quietly = TRUE)) {
    stop("MultiAssayExperiment package required, please install it.",
         call. = FALSE)
  }
  stop("MultiAssayExpermient support not yet implemented")
}
