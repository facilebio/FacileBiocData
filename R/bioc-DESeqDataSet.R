#' @include api.R
NULL

#' @export
#' @noRd
#' @importClassesFrom DESeq2 DESeqDataSet
setClass("FacileDESeqDataSet",
         contains = c("FacileBiocDataStore", "DESeqDataSet"))

#' @rdname FacileBiocDataStore-class
#' @section DESeqDataSet:
#' When a DESeqDataSet object is facilitate()'d, a normalized count matrix
#' will be calculated using `DESeq2::counts(x, normalized = TRUE)` and stored as
#' a matrix named `"normcounts"` in its `assays()` list. These are the values
#' that are returned by (fetch|with)_assay_data when `normalized = TRUE`, which
#' differs from the edgeR:cpm normalized count data which is usually returned
#' from most every other expression container.
#'
#' By default, these normalized counts will be log2 transformed when returned to
#' conform to the expectation in the facilebio ecosystem. To get the deault
#' DESeq2 behaviour, the user would use
#' `fetch_assay_data(.., normalized = TRUE, log = FALSE)`.
#'
#' This function will also look for variance stabilized versions of the
#' data in the `"vst"` and `"rlog"` assay matrices. If no `"vst"` assay is
#' present, it will be run and stored there, unless the `facilitate,run_vst`
#' parameter is set to `FALSE`. This data can be returned using
#' `assay_name = "vst"`
#'
#' @export
#' @param run_vst should we re-run the vst transformation for a DESeqDataSet.
#'   If the `DESeqDataSet` already has a `"vst"` assay, then we'll just take
#'   that, otherwise if this isn't set to `FALSE` it will be run.
#' @param blind,nsub,fitType parameters to send to [DESeq2::vst()] to tweak
#'   how it is run internally
#' @examples
#' # DESeq2 --------------------------------------------------------------------
#' dds <- DESeq2::makeExampleDESeqDataSet(n=2000, m=20)
#' fd <- facilitate(dds)
#' fetch_assay_data(samples(fd), c("gene1", "gene20"))
#' fetch_assay_data(samples(fd), c("gene1", "gene20"), normalized = TRUE)
#'
#' samples(fd) |>
#'   with_assay_data(c("gene1", "gene20"), normalized = TRUE)
#'
#' # Retrieiving different flavors of normalized expression data
#' dat <- samples(fd) |>
#'   with_assay_data("gene1", normalized = TRUE) |>
#'   with_assay_data("gene1", assay_name = "vst") |>
#'   select(-(1:2))
#' colnames(dat) <- c("normcounts", "vst")
#' pairs(dat)
#'
#' \dontrun{
#' dpca <- FacileAnalysis::fpca(fd, assay_name = "vst")
#' FacileAnalysis::shine(dpca)
#' }
facilitate.DESeqDataSet <- function(x, assay_type = "rnaseq",
                                    feature_type = "infer",
                                    organism = "unknown", ...,
                                    run_vst = NULL, blind = TRUE,
                                    nsub = 1000, fitType = "parametric",
                                    verbose = FALSE) {
  reqpkg("DESeq2")
  reqpkg("S4Vectors")
  sinfo <- .init_pdata(x, ...)
  colnames(x) <- rownames(sinfo)
  sinfo <- S4Vectors::DataFrame(sinfo)

  finfo <- .init_fdata(x, feature_type = feature_type, ...)
  ainfo <- .init_assay_info(x, finfo, assay_type = assay_type, ...)

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

  normcounts <- DESeq2::counts(x, normalized = TRUE)
  x <- SummarizedExperiment::`assay<-`(x, "normcounts", value = normcounts)

  anames <- SummarizedExperiment::assayNames(x)
  if ("vst" %in% anames) {
    if (is.null(run_vst)) run_vst <- FALSE
  } else {
    if (is.null(run_vst)) run_vst <- TRUE
  }

  if (isTRUE(run_vst)) {
    if (verbose) {
      if ("vst" %in% anames) {
        message("Rerrunning VST even though it's already stored in ",
                "the container")
      } else {
        message("Running variance stabilizing transform and storing results ",
                "in the 'vst' assay_name")
      }
    }
    vsd <- DESeq2::vst(x, blind = blind, nsub = nsub, fitType = fitType)
    vst <- SummarizedExperiment::assay(vsd)
    x <- SummarizedExperiment::`assay<-`(x, "vst", value = vst)
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
    counts = list(
      assay_type = assay_type,
      feature_type = ainfo[["feature_type"]]),
    normcounts = list(
      assay_type = "tpm",
      feature_type = ainfo[["feature_type"]]),
    vst = list(
      assay_type = "lognorm",
      feature_type = ainfo[["feature_type"]]),
    rlog = list(
      assay_type = "lognorm",
      feature_type = ainfo[["feature_type"]]))
  out@facile[["default_assay"]] <- "counts"
  out@facile[["assay_sample_info"]] <- .init_assay_sample_info(out)
  out@facile[["organism"]] <- organism
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
    normalized = FALSE, batch = NULL, main = NULL, as.matrix = FALSE, ...,
    log = normalized || !is.null(batch), replaced = FALSE, prior.count = 0.1,
    aggregate = FALSE, aggregate.by= "ewm", verbose = FALSE) {
  assert_string(assay_name)
  if (normalized && assay_name == "counts") {
    assay_name <-  "normcounts"
  }
  if (assay_name == "normcounts") {
    normalized <- TRUE
  }
  if (assay_name %in% c("vst", "rlog")) {
    if (!assay_name %in% assay_names(x)) {
      stop("If you want to use vst or rlog data, store the ", assay_name, " ",
           "transformation as an assay in the DESeqDataSet prior to calling ",
           "facilitate()\n. The vst can also be run by calling ",
           "facilitate(x, run_vst = TRUE). See help for more information.")
    }
    normalized <- !is.null(batch)
    log <- FALSE
  }

  fetch_assay_data.FacileBiocDataStore(
    x, features, samples, assay_name = assay_name, normalized = normalized,
    batch = batch, main = main, as.matrix = as.matrix, log = log, ...,
    aggregate = aggregate, aggregate.by = aggregate.by, verbose = verbose)
}

# bioc data retrieval methods --------------------------------------------------

#' @noRd
fdata.DESeqDataSet <- function(x, assay_name = default_assay(x), ...) {
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
assay_names.DESeqDataSet <- function(x, ...) {
  reqpkg("DESeq2")
  SummarizedExperiment::assayNames(x)
}
