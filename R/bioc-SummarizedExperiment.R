#' @include api.R
NULL

#' @export
#' @importFrom utils packageVersion
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
#' @importClassesFrom S4Vectors DataFrame SimpleList
setClass("FacileSummarizedExperiment",
         contains = c("FacileBiocDataStore", "SummarizedExperiment"))

#' @export
#' @noRd
#' @rdname facilitate
#' @method facilitate SummarizedExperiment
facilitate.SummarizedExperiment <- function(x, ...) {

}

# bioc data retrieval methods --------------------------------------------------

#' @noRd
fdata.SummarizedExperiment <- function(x, ...) {
  if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) {
    stop("SummarizedExperiment package required, please install it.",
         call. = FALSE)
  }
  as.data.frame(SummarizedExperiment::rowData(x))
}

#' @noRd
pdata.SummarizedExperiment <- function(x, ...) {
  if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) {
    stop("SummarizedExperiment package required, please install it.",
         call. = FALSE)
  }
  as.data.frame(SummarizedExperiment::colData(x))
}

#' @noRd
adata.SummarizedExperiment <- function(x, name = NULL, ...) {
  if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) {
    stop("SummarizedExperiment package required, please install it.",
         call. = FALSE)
  }
  if (is.null(name)) {
    name <- SummarizedExperiment::assayNames(x)[1L]
  }
  SummarizedExperiment::assay(x, name)
}
