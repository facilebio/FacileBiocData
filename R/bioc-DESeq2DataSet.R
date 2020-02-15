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
  reqpkg("DESeq2")

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
  out@facile[["default_assay"]] <- "counts"
  out
}

# bioc data retrieval methods --------------------------------------------------

#' @noRd
fdata.DESeqDataSet <- function(x, ...) {
  reqpkg("DESeq2")
  as.data.frame(SummarizedExperiment::rowData(x))
}

#' @noRd
pdata.DESeqDataSet <- function(x, ...) {
  reqpkg("DESeq2")
  as.data.frame(SummarizedExperiment::colData(x))
}

#' @noRd
adata.DESeqDataSet <- function(x, name = default_assay(x), ...) {
  reqpkg("DESeq2")
  SummarizedExperiment::assay(x, name)
}

# facile -----------------------------------------------------------------------

#' @noRd
#' @export
assay_names.DESeqDataSet <- function(x, ...) {
  reqpkg("DESeq2")
  SummarizedExperiment::assayNames(x)
}
