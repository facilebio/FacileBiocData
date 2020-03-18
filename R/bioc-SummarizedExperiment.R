#' @include api.R
NULL

#' @export
setClass("FacileSummarizedExperiment",
         contains = c("FacileBiocDataStore", "SummarizedExperiment"))

#' @export
#' @noRd
facilitate.SummarizedExperiment <- function(x, feature_type = "infer",
                                            assay_type = NULL,
                                            assay_info = NULL, ...) {
  reqpkg("SummarizedExperiment")
  is.ranged <- is(x, "RangedSummarizedExperiment") # airway dataset killed me

  if (is.null(SummarizedExperiment::assayNames(x))) {
    anames. <- paste0("adata", seq(SummarizedExperiment::assays(x)) - 1L)
    anames. <- sub("a0$", "a", anames.)
    x <- SummarizedExperiment::`assayNames<-`(x, value = anames.)
  }

  if (test_string(assay_type) && is.null(assay_info)) {
    assay_info <- list(info = list(assay_type = assay_type))
    names(assay_info) <- SummarizedExperiment::assayNames(x)[1L]
  }

  if (is.null(assay_info)) assay_info <- list()

  sinfo <- .init_pdata(x, ...)
  colnames(x) <- rownames(sinfo)
  sinfo <- S4Vectors::DataFrame(sinfo)

  finfo <- .init_fdata(x, ...)
  rownames(x) <- finfo[["feature_id"]]
  finfo <- S4Vectors::DataFrame(finfo)

  x <- SummarizedExperiment::`colData<-`(x, value = sinfo)
  x <- SummarizedExperiment::`rowData<-`(x, value = finfo)

  out <- new(
    "FacileSummarizedExperiment",
    colData = x@colData,
    # Need this to support RangedSummarizedExperiments, like the airway dataset
    assays = SummarizedExperiment::Assays(as(x@assays, "SimpleList")),
    NAMES = x@NAMES,
    elementMetadata = x@elementMetadata,
    metadata = x@metadata)

  if (is.ranged) {
    # Converting the airway dataset really hurt. rownames aren't getting set
    # and the colData is not being passed through ...
    # https://github.com/facilebio/FacileBiocData/issues/2
    out <- SummarizedExperiment::`colData<-`(out, value = sinfo)
    out <- SummarizedExperiment::`rowData<-`(out, value = finfo)
  }

  if (is.null(rownames(out))) {
    # I still don't get how this happens, since rownames(x) was set
    rownames(out) <- finfo[["feature_id"]]
  }
  out@facile[["assay_info"]] <- assay_info
  out@facile[["default_assay"]] <- SummarizedExperiment::assayNames(out)[1L]
  out@facile[["assay_sample_info"]] <- .init_assay_sample_info(out)
  out
}

# #' @export
# #' @noRd
# facilitate.RangedSummarizedExperiment <- function(x, assay_type = "infer",
#                                                   feature_type = "infer", ...) {
#   reqpkg("SummarizedExperiment")
#   if (is.null(SummarizedExperiment::assayNames(x))) {
#     anames. <- paste0("adata", seq(SummarizedExperiment::assays(x)) - 1L)
#     anames. <- sub("a0$", "a", anames.)
#     x <- SummarizedExperiment::`assayNames<-`(x, value = anames.)
#   }
#
#   sinfo <- .init_pdata(x, ...)
#   colnames(x) <- rownames(sinfo)
#   sinfo <- S4Vectors::DataFrame(sinfo)
#
#   finfo <- .init_fdata(x, ...)
#   rownames(x) <- finfo[["feature_id"]]
#   finfo <- S4Vectors::DataFrame(finfo)
#
#   x <- SummarizedExperiment::`colData<-`(x, value = sinfo)
#   x <- SummarizedExperiment::`rowData<-`(x, value = finfo)
#
#
#   out <- new("FacileSummarizedExperiment",
#              colData = x@colData,
#              assays = x@assays,
#              NAMES = x@NAMES,
#              elementMetadata = x@elementMetadata,
#              metadata = x@metadata)
#
#   out@facile[["assay_sample_info"]] <- .init_assay_sample_info(out)
#   out
# }


# bioc data retrieval methods --------------------------------------------------

#' @noRd
fdata.SummarizedExperiment <- function(x, assay_name = default_assay(x), ...) {
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
  out <- SummarizedExperiment::assay(x, name)
  .cleanup_adata(x, out, name = name, ...)
}

# facile -----------------------------------------------------------------------

#' @noRd
#' @export
assay_names.SummarizedExperiment <- function(x, ...) {
  reqpkg("SummarizedExperiment")
  SummarizedExperiment::assayNames(x)
}
