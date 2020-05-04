#' @include api.R
#' @include bioc-DGEList.R
NULL

#' @export
setClass("FacileEList", contains = c("FacileBiocDataStore", "EList"))

#' @export
#' @noRd
facilitate.EList <- function(x, assay_type = "lognorm", feature_type = "infer",
                             organism = "unknown", ...) {
  reqpkg("limma")

  sinfo <- .init_pdata(x, ...)
  colnames(x) <- rownames(sinfo)

  E <- x[["E"]]
  colnames(E) <- rownames(sinfo)

  # Currently we only support one assay
  finfo <- .init_fdata(x, feature_type = feature_type, ...)
  rownames(E) <- finfo[["feature_id"]]

  x$E <- E
  x$targets <- sinfo
  x$genes <- finfo

  ainfo <- .init_assay_info(x, finfo, assay_type = assay_type, ...)

  out <- new("FacileEList", lapply(x, identity))
  out@facile[["assay_info"]] <- list(
    E = list(assay_type = assay_type, feature_type = ainfo[["feature_type"]]))
  out@facile[["default_assay"]] <- "E"
  out@facile[["assay_sample_info"]] <- .init_assay_sample_info(out)
  out@facile[["organism"]] <- organism
  out
}

# bioc data retrieval methods --------------------------------------------------

#' @noRd
fdata.EList <- function(x, assay_name = default_assay(x), ...) {
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
