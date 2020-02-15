#' @include api.R
NULL

#' @export
setClass("FacileDGEList", contains = c("FacileBiocDataStore", "DGEList"))

#' @export
#' @noRd
#' @rdname facilitate
#' @method facilitate DGEList
facilitate.DGEList <- function(x, ...) {
  reqpkg("edgeR")
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
  reqpkg("edgeR")
  x[["genes"]]
}

#' @noRd
pdata.DGEList <- function(x, ...) {
  reqpkg("edgeR")
  x[["samples"]]
}

#' @noRd
adata.DGEList <- function(x, name = default_assay(x), ...) {
  reqpkg("edgeR")
  if (is.null(name)) {
    name <- "counts"
  }
  out <- x[[name]]
  assert_matrix(out, "numeric", nrows = nrow(x), ncols = ncol(x))
  out
}

# facile -----------------------------------------------------------------------

#' @noRd
#' @export
assay_names.DGEList <- function(x, ..., .required = "counts") {
  reqpkg("edgeR")
  dnames <- names(x)
  if (!.required %in% dnames) {
    stop("`counts` not found in DGEList, you've got problems")
  }
  out <- unique(c(.required, dnames))
  for (dname in dnames) {
    elem <- x[[dname]]
    if (!(is.matrix(elem) && nrow(elem) == nrow(x) && ncol(elem) == ncol(x))) {
      out <- setdiff(out, dname)
    }
  }
  out
}
