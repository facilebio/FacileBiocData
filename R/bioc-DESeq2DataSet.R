#' @include api.R
NULL

#' @export
#' @importClassesFrom DESeq2 DESeqDataSet
setClass("FacileDESeqDataSet",
         contains = c("FacileBiocDataStore", "DESeqDataSet"))

#' @export
#' @noRd
#' @rdname facilitate
#' @method facilitate SummarizedExperiment
facilitate.DESeqDataSet <- function(x, assay_name = NULL, ...) {
  if (!requireNamespace("DESeq2", quietly = TRUE)) {
    stop("DESeq2 package required, please install it.",
         call. = FALSE)
  }

  sinfo <- .init_pdata(x, ...)
  colnames(x) <- sinfo[["sample_id"]]
  sinfo <- S4Vectors::DataFrame(sinfo)

  finfo <- .init_fdata(x, assay_name = assay_name, ...)
  # because of dds@rowRanges, we can't have GRanges-like names in here
  axe.cols <- c("seqnames", "ranges", "strand", "start", "end", "width",
                "element")
  axe.cols <- intersect(colnames(finfo), axe.cols)
  if (length(axe.cols)) {
    finfo <- finfo[, !colnames(finfo) %in% axe.cols]
  }
  rownames(x) <- finfo[["feature_id"]]
  finfo <- S4Vectors::DataFrame(finfo)

  x <- SummarizedExperiment::`colData<-`(x, value = sinfo)
  x <- SummarizedExperiment::`rowData<-`(x, value = finfo)


  out <- new("FacileDESeqDataSet",
             design = x@design,
             dispersionFunction = x@dispersionFunction,
             rowRanges = x@rowRanges,
             colData = x@colData,
             assays = x@assays,
             NAMES = x@NAMES,
             elementMetadata = x@elementMetadata,
             metadata = x@metadata)
  out
}

# bioc data retrieval methods --------------------------------------------------

#' @noRd
fdata.DESeqDataSet <- function(x, ...) {
  if (!requireNamespace("DESeq2", quietly = TRUE)) {
    stop("DESeq2 package required, please install it.",
         call. = FALSE)
  }
  as.data.frame(SummarizedExperiment::rowData(x))
}

#' @noRd
pdata.DESeqDataSet <- function(x, ...) {
  if (!requireNamespace("DESeq2", quietly = TRUE)) {
    stop("SummarizedExperiment package required, please install it.",
         call. = FALSE)
  }
  as.data.frame(SummarizedExperiment::colData(x))
}

#' @noRd
adata.DESeqDataSet <- function(x, name = NULL, ...) {
  if (!requireNamespace("DESeq2", quietly = TRUE)) {
    stop("SummarizedExperiment package required, please install it.",
         call. = FALSE)
  }
  if (is.null(name)) {
    name <- SummarizedExperiment::assayNames(x)[1L]
  }
  SummarizedExperiment::assay(x, name)
}

#' @noRd
anames.DESeqDataSet <- function(x, ...) {
  if (!requireNamespace("DESeq2", quietly = TRUE)) {
    stop("SummarizedExperiment package required, please install it.",
         call. = FALSE)
  }
  SummarizedExperiment::assayNames(x)
}
