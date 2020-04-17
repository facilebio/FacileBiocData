#' @include api.R
NULL

#' @export
#' @noRd
#' @importClassesFrom DESeq2 DESeqDataSet
setClass("FacileDESeqDataSet",
         contains = c("FacileBiocDataStore", "DESeqDataSet"))

#' @rdname FacileBiocDataStore
#' @section DESeqDataSet:
#'
#' The FacileDESeqDataSet will look for variance stabilized versions of the
#' data in the `"vst"` and `"rlog"` assay matrices. If no `"vst"` assay is
#' present, it will be run and stored there, unless the `facilitate,run_vst`
#' parameter is set to `FALSE`.
#'
#' Because DESeq uses a different normalization method than edgeR's TMM, when
#' the user calls `fetch_assay_data(.., normalized = TRUE)`, the default will
#' be to return the normalized count data retrieved from
#' [DESeq2::counts()] with `normalized = TRUE`.
#'
#' To return [edgeR::cpm()] values, you can set `normalized = "cpm"`, but this
#' must be working over the `"counts"` assay.
#'
#' 1. Add parameters to run vst/rlog?
#' 2. Enable vst, rlog, and normcounts  to be retrieved via
#'    fetch_assay_data(assay_name = {"vst"|"rlog"|"normcounts"})
#'
#' @export
#' @examples
#'
#' # DESeq2 --------------------------------------------------------------------
#' dds <- DESeq2::makeExampleDESeqDataSet(n=2000, m=20)
#' fd <- facilitate(dds)
#' fetch_assay_data(samples(fd), c("gene1", "gene20"))
#' fetch_assay_data(samples(fd), c("gene1", "gene20"), normalized = TRUE)
#' fetch_assay_data(samples(fd), c("gene1", "gene20"), normalized = "cpm")
#'
#' samples(fd) %>%
#'   with_assay_data(c("gene1", "gene20"), normalized = TRUE)
#'
#' # Retrieiving different flavors of normalized expression data
#' dat <- samples(fd) %>%
#'   with_assay_data("gene1", normalized = TRUE) %>%
#'   with_assay_data("gene1", normalized = "cpm") %>%
#'   with_assay_data("gene1", assay_name = "vst") %>%
#'   select(-(1:2))
#' colnames(dat) <- c("normcounts", "cpm", "vst")
#' pairs(dat)
#'
#' dpca <- FacileAnalysis::fpca(fd, assay_name = "vst")
facilitate.DESeqDataSet <- function(x, assay_type = "rnaseq",
                                    feature_type = "infer", ...,
                                    run_vst = NULL, blind = TRUE,
                                    nsub = 1000, fitType = "parametric",
                                    prior.count = 0.1, verbose = FALSE) {
  reqpkg("DESeq2")

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
    vst = list(
      assay_type = "lognorm",
      feature_type = ainfo[["feature_type"]]),
    rlog = list(
      assay_type = "lognorm",
      feature_type = ainfo[["feature_type"]]),
    normcounts = list(
      assay_type = "lognorm",
      feature_type = ainfo[["feature_type"]]))
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
  if (test_string(normalized) && normalized == "cpm") {
    # This is for DESeqDataSet that wans to use edgeR normalized counts
    normalized <- TRUE
    assay_name <- "counts"
  } else if (assay_name == "counts") {
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

  fetch_assay_data.FacileBiocDataStore(
    x, features, samples, assay_name = assay_name, normalized = normalized,
    batch = batch, main = main, as.matrix = as.matrix, ...,
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
