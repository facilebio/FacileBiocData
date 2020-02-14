#' @include api.R
NULL

# NOTE: FacileSingleCellExperiment should be a subclass of FacileSumarrizedExperiment,
#       figure it out!
#' @export
#' @importFrom utils packageVersion
#' @importClassesFrom SingleCellExperiment SingleCellExperiment
#' @importClassesFrom S4Vectors DataFrame SimpleList
setClass("FacileSingleCellExperiment",
         contains = c("FacileBiocDataStore", "SingleCellExperiment"))

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
