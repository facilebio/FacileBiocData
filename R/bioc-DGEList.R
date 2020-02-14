#' @include api.R
NULL

#' @export
setClass("FacileDGEList",
         slots = c(facile = "list"),
         contains = c("FacileBiocDataStore", "DGEList"),
         prototype = prototype(facile = list()))

#' @export
#' @noRd
#' @rdname facilitate
#' @method facilitate DGEList
facilitate.DGEList <- function(x, ...) {
  sinfo <- .init_pdata(x, ...)
  counts <- x[["counts"]]
  colnames(counts) <- rownames(sinfo)
  out <- new("FacileDGEList",
             list(counts = counts,
                  samples = sinfo,
                  genes = x[["genes"]]))
  # eav <- as.EAVtable(sinfo)
  # out@facile[["eav"]] <- eav
  # out@facile[["covariate_def"]] <- attr(eav, "covariate_def")
  out
}

# bioc data retrieval methods --------------------------------------------------

#' @noRd
fdata.DGEList <- function(x, ...) {
  x[["genes"]]
}

#' @noRd
pdata.DGEList <- function(x, ...) {
  x[["samples"]]
}

#' @noRd
adata.DGEList <- function(x, name = "counts", ...) {
  out <- x[[name]]
  assert_matrix(out, "numeric", nrows = nrow(fdata(x)), ncols = ncol(pdata(x)))
  out
}
