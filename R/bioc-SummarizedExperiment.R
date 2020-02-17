#' @include api.R
NULL

#' @export
setClass("FacileSummarizedExperiment",
         contains = c("FacileBiocDataStore", "SummarizedExperiment"))

#' @export
#' @noRd
#' @rdname facilitate
#' @method facilitate SummarizedExperiment
facilitate.SummarizedExperiment <- function(x, assay_type = "infer",
                                            feature_type = "infer", ...) {
  reqpkg("SummarizedExperiment")
  if (is.null(SummarizedExperiment::assayNames(x))) {
    anames. <- paste0("adata", seq(SummarizedExperiment::assays(x)) - 1L)
    anames. <- sub("a0$", "a", anames.)
    x <- SummarizedExperiment::`assayNames<-`(x, value = anames.)
  }

  sinfo <- .init_pdata(x, ...)
  colnames(x) <- sinfo[["sample_id"]]
  sinfo <- S4Vectors::DataFrame(sinfo)

  finfo <- .init_fdata(x, ...)
  rownames(x) <- finfo[["feature_id"]]
  finfo <- S4Vectors::DataFrame(finfo)

  x <- SummarizedExperiment::`colData<-`(x, value = sinfo)
  x <- SummarizedExperiment::`rowData<-`(x, value = finfo)


  out <- new("FacileSummarizedExperiment",
             colData = x@colData,
             assays = x@assays,
             NAMES = x@NAMES,
             elementMetadata = x@elementMetadata,
             metadata = x@metadata)

  out@facile[["assay_sample_info"]] <- .init_assay_sample_info(out)
  out
}

# bioc data retrieval methods --------------------------------------------------

#' @noRd
fdata.SummarizedExperiment <- function(x, ...) {
  reqpkg("SummarizedExperiment")
  as.data.frame(SummarizedExperiment::rowData(x))
}

#' @noRd
pdata.SummarizedExperiment <- function(x, ...) {
  reqpkg("SummarizedExperiment")
  as.data.frame(SummarizedExperiment::colData(x))
}

#' @noRd
adata.SummarizedExperiment <- function(x, name = default_assay(x), ...) {
  reqpkg("SummarizedExperiment")
  if (is.null(name)) {
    name <- SummarizedExperiment::assayNames(x)[1L]
  }
  SummarizedExperiment::assay(x, name)
}

# facile -----------------------------------------------------------------------

#' @noRd
#' @export
assay_names.SummarizedExperiment <- function(x, ...) {
  reqpkg("SummarizedExperiment")
  SummarizedExperiment::assayNames(x)
}
