#' @include api.R
NULL

#' @export
#' @noRd
#' @method facilitate DGEList
facilitate.EList <- function(x, ...) {

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
  x[["genes"]]
}

#' @noRd
pdata.EList <- function(x, ...) {
  x[["targets"]]
}

#' @noRd
adata.EList <- function(x, name = "E", ...) {
  out <- x[[name]]
  assert_matrix(out, "numeric", nrows = nrow(fdata(x)), ncols = ncol(pdata(x)))
  out
}
