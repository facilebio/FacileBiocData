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
  reqpkg("SingleCellExperiment")
  stop("SingleCellExperiment support not yet implemented")
}

# bioc data retrieval methods --------------------------------------------------

#' @noRd
fdata.SingleCellExperiment <- function(x, ...) {
  reqpkg("SingleCellExperiment")
  as.data.frame(SingleCellExperiment::rowData(x))
}

#' @noRd
pdata.SingleCellExperiment <- function(x, ...) {
  reqpkg("SingleCellExperiment")
  as.data.frame(SingleCellExperiment::colData(x))
}

#' @noRd
adata.SingleCellExperiment <- function(x, name = default_assay(x), ...) {
  reqpkg("SingleCellExperiment")
  SingleCellExperiment::assay(x, name)
}

# facile -----------------------------------------------------------------------

#' @noRd
#' @export
assay_names.SingleCellExperiment <- function(x, ...) {
  reqpkg("SingleCellExperiment")
  stop("SingleCellExperiment support not yet implemented")
}

