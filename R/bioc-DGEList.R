#' @include api.R
NULL

#' @export
setClass("FacileDGEList", contains = c("FacileBiocDataStore", "DGEList"))

#' @rdname FacileBiocDataStore
#' @section DGEList:
#' We assume the DGEList holds `"rnaseq"` assay data. Set the `assay_type`
#' parameter if that's not the case.
#'
#' @export
#' @examples
#' # edgeR ---------------------------------------------------------------------
#' y <- example_bioc_data(class = "DGEList")
#' yf <- facilitate(y)
#' fpca(yf)
facilitate.DGEList <- function(x, assay_type = "rnaseq", ...) {
  reqpkg("edgeR")
  sinfo <- .init_pdata(x, ...)
  colnames(x) <- rownames(sinfo)

  counts <- x[["counts"]]
  colnames(counts) <- rownames(sinfo)

  # Currently we only support one assay
  finfo <- .init_fdata(x, assay_type = assay_type, ...)
  rownames(counts) <- finfo[["feature_id"]]

  x$counts <- counts
  x$samples <- sinfo
  x$genes <- finfo
  out <- new("FacileDGEList", lapply(x, identity))
  out@facile[["assay_info"]] <- list(counts = list(assay_type = assay_type))
  out@facile[["default_assay"]] <- "counts"
  out@facile[["assay_sample_info"]] <- .init_assay_sample_info(out)
  out
}

# bioc data retrieval methods --------------------------------------------------

#' @noRd
fdata.DGEList <- function(x, assay_name = default_assay(x), ...) {
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
  .cleanup_adata(x, out, name = name, ...)
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
