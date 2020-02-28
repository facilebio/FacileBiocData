#' @include api.R
NULL

# Why don't we have to `@importClassesFrom Biobase ExpressionSet` here like we
# do for DESeqDataSet?

#' @export
setClass("FacileExpressionSet",
         contains = c("FacileBiocDataStore", "ExpressionSet"))

#' @export
#' @noRd
facilitate.ExpressionSet <- function(x, assay_type = "infer",
                                     feature_type = "infer", ...) {
  if (!requireNamespace("Biobase", quietly = TRUE)) {
    stop("Biobase package required, please install it.",
         call. = FALSE)
  }

  out <- new("FacileExpressionSet",
             experimentData = x@experimentData,
             assayData = x@assayData,
             phenoData = x@phenoData,
             featureData = x@featureData,
             anotation = x@annotation,
             protocolData = x@protocolData)

  sinfo <- .init_pdata(x, ...)
  colnames(x) <- rownames(sinfo)
  out <- Biobase::`pData<-`(out, sinfo)

  # Currently we only support one assay
  finfo <- .init_fdata(x, ...)
  rownames(out) <- finfo[["feature_id"]]
  out <- Biobase::`fData<-`(out, finfo)
  # eav <- as.EAVtable(sinfo)
  # out@facile[["eav"]] <- eav
  # out@facile[["covariate_def"]] <- attr(eav, "covariate_def")
  out@facile[["assay_sample_info"]] <- .init_assay_sample_info(out)
  out
}

#' @noRd
ifacile.FacileExpressionSet <- function(x, ...) x@facile

# bioc data retrieval methods --------------------------------------------------

#' @noRd
fdata.ExpressionSet <- function(x, ...) {
  if (!requireNamespace("Biobase", quietly = TRUE)) {
    stop("Biobase package required, please install it.",
         call. = FALSE)
  }
  Biobase::fData(x)
}

#' @noRd
pdata.ExpressionSet <- function(x, ...) {
  reqpkg("Biobase")
  Biobase::pData(x)
}

#' @noRd
adata.ExpressionSet <- function(x, name = default_assay(x), ...) {
  reqpkg("Biobase")
  if (is.null(name)) {
    name <- Biobase::assayDataElementNames(eset)[1L]
  }
  Biobase::assayDataElement(x, name)
}

# facile -----------------------------------------------------------------------

#' @noRd
#' @export
assay_names.ExpressionSet <- function(x, ...) {
  reqpkg("Biobase")
  out <- Biobase::assayDataElementNames(x)
  if ("exprs" %in% out) {
    out <- c("exprs", setdiff(out, "exprs"))
  }
  out
}
