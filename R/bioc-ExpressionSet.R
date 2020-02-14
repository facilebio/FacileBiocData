#' @include api.R
NULL

# Why don't we have to `@importClassesFrom Biobase ExpressionSet` here like we
# do for DESeqDataSet?

#' @export
setClass("FacileExpressionSet",
         contains = c("FacileBiocDataStore", "ExpressionSet"))

#' @export
#' @noRd
#' @rdname facilitate
#' @method facilitate ExpressionSet
facilitate.ExpressionSet <- function(x, assay_name = NULL, ...) {
  if (!requireNamespace("Biobase", quietly = TRUE)) {
    stop("Biobase package required, please install it.",
         call. = FALSE)
  }

  # out <- new("FacileExpressionSet",
  #            experimentData = Biobase::experimentData(x),
  #            assayData = Biobase::assayData(x),
  #            phenoData = Biobase::phenoData(x),
  #            featureData = Biobase::featureData(x),
  #            anotation = Biobase::annotation(x),
  #            protocolData = Biobase::protocolData(x))
  #
  # We are already digging directly into the slotNames of the ExpressionSet
  # class, so just use the `@` accessor directly
  out <- new("FacileExpressionSet",
             experimentData = x@experimentData,
             assayData = x@assayData,
             phenoData = x@phenoData,
             featureData = x@featureData,
             anotation = x@annotation,
             protocolData = x@protocolData)

  sinfo <- .init_pdata(x, ...)
  colnames(out) <- sinfo[["sample_id"]]
  out <- Biobase::`pData<-`(out, sinfo)

  # Currently we only support one assay
  finfo <- .init_fdata(x, assay_name = assay_name, ...)
  rownames(out) <- finfo[["feature_id"]]
  out <- Biobase::`fData<-`(out, finfo)
  # eav <- as.EAVtable(sinfo)
  # out@facile[["eav"]] <- eav
  # out@facile[["covariate_def"]] <- attr(eav, "covariate_def")
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
  if (!requireNamespace("Biobase", quietly = TRUE)) {
    stop("Biobase package required, please install it.",
         call. = FALSE)
  }
  Biobase::pData(x)
}

#' @noRd
adata.ExpressionSet <- function(x, name = "exprs", ...) {
  if (!requireNamespace("Biobase", quietly = TRUE)) {
    stop("Biobase package required, please install it.",
         call. = FALSE)
  }
  if (is.null(name)) {
    name <- Biobase::assayDataElementNames(eset)[1L]
  }
  Biobase::assayDataElementy(x, name)
}

#' @noRd
anames.ExpressionSet <- function(x, ...) {
  if (!requireNamespace("Biobase", quietly = TRUE)) {
    stop("Biobase package required, please install it.",
         call. = FALSE)
  }
  Biobase::assayDataElementNames(x)
}
