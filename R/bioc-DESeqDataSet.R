#' @include api.R
NULL

#' @export
#' @importClassesFrom DESeq2 DESeqDataSet
setClass("FacileDESeqDataSet",
         contains = c("FacileBiocDataStore", "DESeqDataSet"))

#' @section DESeqDataSet
#' 1. Add parameters to run vst/rlog?
#' 2. Enable vst, rlog, and normcounts  to be retrieved via
#'    fetch_assay_data(assay_name = {"vst"|"rlog"|"normcounts"})
#'
#' @export
#' @noRd
#' @examples
#' dds <- DESeq2::makeExampleDESeqDataSet(n=2000, m=20)
#' fd <- facilitate(dds)
#' fetch_assay_data(samples(fd), c("gene1", "gene20"))
#' fetch_assay_data(samples(fd), c("gene1", "gene20"), normalized = TRUE)
#' fetch_assay_data(samples(fd), c("gene1", "gene20"), normalized = TRUE,
#'                  assay_name = "cpm")
#' samples(fd) %>%
#'   with_assay_data(c("gene1", "gene20"), normalized = TRUE)
#'
#' if (FALSE) {
#' dat <- samples(fd) %>%
#'   with_assay_data("gene1", normalized = TRUE) %>%
#'   with_assay_data("gene1", assay_name = "vst") %>%
#'   with_assay_data("gene1", assay_name = "cpm")
#' pairs(as.matrix(dat[, -(1:2)]))
#' abline(0, 1, col = "red")
#' }
#'
#'
facilitate.DESeqDataSet <- function(x, assay_type = "rnaseq",
                                    feature_type = "infer", ...,
                                    run_vst = TRUE, blind = TRUE,
                                    nsub = 1000, fitType = "parametric") {
  reqpkg("DESeq2")

  sinfo <- .init_pdata(x, ...)
  colnames(x) <- rownames(sinfo)
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

  if (is.null(DESeq2::sizeFactors(x))) {
    x <- DESeq2::estimateSizeFactors(x)
  }

  anames <- SummarizedExperiment::assayNames(x)
  if (!"vst" %in% anames) {
    if (missing(run_vst)) {
      message("Running vst on dataset, next time either include a vst assay, ",
              "or set `run_vst = FALSE`")
    }
    if (run_vst) {
      vsd <- DESeq2::vst(x, blind = blind, nsub = nsub, fitType = fitType)
      vst <- SummarizedExperiment::assay(vsd)
      x <- SummarizedExperiment::`assay<-`(x, "vst", value = vst)
    }
  }

  out <- new("FacileDESeqDataSet",
             design = x@design,
             dispersionFunction = x@dispersionFunction,
             rowRanges = x@rowRanges,
             colData = x@colData,
             assays = x@assays,
             NAMES = x@NAMES,
             elementMetadata = x@elementMetadata,
             metadata = x@metadata)
  out@facile[["assay_info"]] <- list(
    counts = list(assay_type = assay_type),
    vst = list(assay_type = "lognorm"),
    rlog = list(assay_type = "lognorm"),
    normcounts = list(assay_type = "lognorm"))
  out@facile[["default_assay"]] <- "counts"
  out@facile[["assay_sample_info"]] <- .init_assay_sample_info(out)
  out
}


#' DESeqDataSet specific fetch_assay_data functions enables us to provide
#' DESeq2 style normalized counts (via `log2(counts(x, normalized = TRUE))`
#' when `normalized = TRUE`, or to intercept the special `"vst"` and `"rlog"`
#' assay_names to use those transformed data for analysis.
#'
#' To get edgeR::cpm style normalized counts, use `assay_name = "cpm"`
#'
#' @noRd
#' @export
fetch_assay_data.FacileDESeqDataSet <- function(
    x, features = NULL, samples = NULL, assay_name = default_assay(x),
    normalized = FALSE, batch = NULL, main = NULL, as.matrix = FALSE,
    ..., prior.count = 0.1, aggregate = FALSE, aggregate.by= "ewm",
    verbose = FALSE) {
  assert_string(assay_name)
  if (assay_name == "counts") {
    if (normalized) {
      assert_number(prior.count, lower = 0.001)
      nc <- log2(DESeq2::counts(x, normalized = TRUE) + prior.count)
      x <- SummarizedExperiment::`assay<-`(x, "normcounts", value = nc)
      assay_name <- "normcounts"
    }
  }
  if (assay_name %in% c("vst", "rlog")) {
    if (!assay_name %in% assay_names(x)) {
      stop("If you want to use vst or rlog data, store the ", assay_name, " ",
           "transformation as an assay in the DESeqDataSet prior to calling ",
           "facilitate()\n. The vst can also be run by calling ",
           "facilitate(x, run_vst = TRUE). See help for more information.")
    }
    normalized <- TRUE
  }
  if (assay_name == "cpm") {
    normalized <- TRUE
    assay_name <- "counts"
  }

  fetch_assay_data.FacileBiocDataStore(
    x, features, samples, assay_name = assay_name, normalized = normalized,
    batch = batch, main = main, as.matrix = as.matrix, ...,
    aggregate = aggregate, aggregate.by = aggregate.by, verbose = verbose)
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
  out <- SummarizedExperiment::assay(x, name)
  .cleanup_adata(x, out, name = name, ...)
}

# facile -----------------------------------------------------------------------

#' @noRd
#' @export
assay_names.DESeqDataSet <- function(x, ..., internal. = FALSE) {
  reqpkg("DESeq2")
  # unique(c(SummarizedExperiment::assayNames(x), c("cpm")))
  if (internal.) {
    SummarizedExperiment::assayNames(x)
  } else {
    unique(c(SummarizedExperiment::assayNames(x), c("cpm")))
  }
}
