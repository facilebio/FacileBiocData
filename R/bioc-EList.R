#' @include api.R
#' @include bioc-DGEList.R
NULL

#' @export
setClass("FacileEList", contains = c("FacileBiocDataStore", "EList"))

#' @export
#' @noRd
facilitate.EList <- function(x, feature_type = "infer", assay_type = "lognorm",
                             ...) {
  reqpkg("limma")

  sinfo <- .init_pdata(x, ...)
  colnames(x) <- rownames(sinfo)

  E <- x[["E"]]
  colnames(E) <- rownames(sinfo)

  # Currently we only support one assay
  finfo <- .init_fdata(x, ...)
  rownames(E) <- finfo[["feature_id"]]

  x$E <- E
  x$targets <- sinfo
  x$genes <- finfo

  out <- new("FacileEList", lapply(x, identity))
  out@facile[["assay_info"]] <- list(E = list(assay_type = assay_type))
  out@facile[["assay_sample_info"]] <- .init_assay_sample_info(out)
  out
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
  assert_matrix(out, "numeric", nrows = nrow(x), ncols = ncol(x))
  .cleanup_adata(x, out, name = name, ...)
}

# facile -----------------------------------------------------------------------

#' @noRd
#' @export
assay_names.EList <- function(x, ...) {
  reqpkg("limma")
  assay_names.DGEList(x, ..., .required = "E")
}