#' @include api.R
NULL

#' @export
#' @importClassesFrom DESeq2 DESeqDataSet
setClass("FacileDESeqDataSet",
         contains = c("FacileBiocDataStore", "DESeqDataSet"))

#' @export
#' @noRd
facilitate.DESeqDataSet <- function(x, assay_type = "rnaseq",
                                    feature_type = "infer", ...) {
  reqpkg("DESeq2")

  sinfo <- .init_pdata(x, ...)
  colnames(x) <- sinfo[["sample_id"]]
  sinfo <- S4Vectors::DataFrame(sinfo)

  finfo <- .init_fdata(x, ...)
  # This was first developed during bioc3.6, and at the time a DESeqDataSet had
  # a GRanges-like @rowRanges slot, which prohibits granges-like colnames in our
  # feature information.
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
  out@facile[["assay_sample_info"]] <- .init_assay_sample_info(out)
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
