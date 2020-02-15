#' @include api.R
#' @include bioc-DGEList.R
NULL

#' @export
setClass("FacileEList", contains = c("FacileBiocDataStore", "EList"))

#' @export
#' @noRd
#' @method facilitate DGEList
facilitate.EList <- function(x, ...) {
  reqpkg("limma")
  stop("facilitate.EList not yet implemented")
  x$facile <- list(
    extra = NULL,
    stuff = NULL,
    here = NULL)

  class(x) <- c("FacileEList", "FacileBiocDataStore", "FacileBiocDataStore",
                class(x))
  x
}

# bioc data retrieval methods --------------------------------------------------

#' @noRd
fdata.EList <- function(x, ...) {
  reqpkg("limma")
  x[["genes"]]
}

#' @noRd
pdata.EList <- function(x, ...) {
  reqpkg("limma")
  x[["targets"]]
}

#' @noRd
adata.EList <- function(x, name = default_assay(x), ...) {
  reqpkg("limma")
  out <- x[[name]]
  assert_matrix(out, "numeric", nrows = nrow(fdata(x)), ncols = ncol(pdata(x)))
  out
}

# facile -----------------------------------------------------------------------

#' @noRd
#' @export
assay_names.EList <- function(x, ...) {
  reqpkg("limma")
  anames.DGEList(x, ..., .required = "E")
}
