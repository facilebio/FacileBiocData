#' @include api.R
NULL

# NOTE: FacileSingleCellExperiment should be a subclass of FacileSumarrizedExperiment,
#       figure it out!
#' @export
#' @importFrom utils packageVersion
#' @importClassesFrom SingleCellExperiment SingleCellExperiment
setClass("FacileSingleCellExperiment",
         contains = c("FacileBiocDataStore", "SingleCellExperiment"))

#' @export
#' @noRd
facilitate.SingleCellExperiment <- function(x, ...) {
  reqpkg("SingleCellExperiment")
  stop("SingleCellExperiment support not yet implemented")
  if (is.null(assay_info)) {
    assay_info <- list(counts = list(assay_type = "rnaseq_sparse"),
                       logcounts = list(assay_type = "lognorm"))
  }
  # I thing "logcounts" would be the default assay instead of "counts" because
  # by the time someone would want to "facilitate" an SCE, it would have already
  # gone through the pre-processing steps.
  out@facile[["default_assay"]] <- "logcounts"
}

# bioc data retrieval methods --------------------------------------------------

#' @noRd
fdata.SingleCellExperiment <- function(x, assay_name = default_assay(x), ...) {
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
  out <- SingleCellExperiment::assay(x, name)
  .cleanup_adata(x, out, name = name, ...)
}

# facile -----------------------------------------------------------------------

#' @noRd
#' @export
assay_names.SingleCellExperiment <- function(x, ...) {
  reqpkg("SingleCellExperiment")
  stop("SingleCellExperiment support not yet implemented")
}

