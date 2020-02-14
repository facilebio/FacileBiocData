#' @include api.R
NULL

#' @export
#' @importFrom utils packageVersion
#' @importClassesFrom SingleCellExperiment SingleCellExperiment
#' @importClassesFrom S4Vectors DataFrame SimpleList
setClass("FacileSingleCellExperiment",
         slots = c(facile = "list"),
         contains = c("FacileBiocDataStore", "SingleCellExperiment"),
         prototype = prototype(facile = list()))

#' @export
#' @noRd
#' @rdname facilitate
#' @method facilitate SingleCellExperiment
facilitate.SingleCellExperiment <- function(x, ...) {

}

# bioc data retrieval methods --------------------------------------------------

#' @noRd
fdata.SingleCellExperiment <- function(x, ...) {
  if (!requireNamespace("SingleCellExperiment", quietly = TRUE)) {
    stop("SingleCellExperiment package required, please install it.",
         call. = FALSE)
  }
  as.data.frame(SingleCellExperiment::rowData(x))
}

#' @noRd
pdata.SingleCellExperiment <- function(x, ...) {
  if (!requireNamespace("SingleCellExperiment", quietly = TRUE)) {
    stop("SingleCellExperiment package required, please install it.",
         call. = FALSE)
  }
  as.data.frame(SingleCellExperiment::colData(x))
}

#' @noRd
adata.SingleCellExperiment <- function(x, name = NULL, ...) {
  if (!requireNamespace("SingleCellExperiment", quietly = TRUE)) {
    stop("SingleCellExperiment package required, please install it.",
         call. = FALSE)
  }
  SingleCellExperiment::assay(x, name)
}
