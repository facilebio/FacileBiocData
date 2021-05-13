#' Map the Bioconductor assay container class to the package it comes from.
#' @noRd
.bioc_types <- tribble(
  ~class,                 ~package,
  "DGEList",              "edgeR",
  "EList",                "limma",
  "ExpressionSet",        "Biobase",
  "SummarizedExperiment", "SummarizedExperiment",
  "DESeqDataSet",         "DESeq2",
  "SingleCellExperiment", "SingleCellExperiment",
  "MultiAssayExperiment", "MultiAssayExperiment")


#' Create an example Bioconductor assay container
#'
#' We can construct the following types of containers:
#'
#' 1. DGEList
#' 2. EList
#' 3. ExpressionSet
#' 4. SummarizedExperiment
#' 5. DESeqDataSet
#' 6. SingleCellExperiment (not yet)
#' 7. MultiAssayExperiment (not yet)
#'
#' @export
#' @param type The type of assay container. Pick from the ones enumerated above.
#' @param efds An already loaded example FacileDataStore to convert to a
#'   Bioconductor container.
#' @param y A DGEList for direct conversion to another Bioconductor object.
#'   This does not to anything with any FacileDataStore anywhere.
#' @return The Bioconductor container
example_bioc_data <- function(class = "DGEList", efds = NULL, y = NULL, ...) {
  info <- filter(.bioc_types, .data$class == .env$class)
  if (nrow(info) != 1L) stop("Unsupported data class:", class, call. = FALSE)
  reqpkg(info[["package"]])

  if (!is(y, "DGEList")) {
    if (!is(efds, "FacileDataSet")) efds <- FacileData::exampleFacileDataSet()
    y <- biocbox(efds, "DGEList")
    y[["samples"]][["samid"]] <- NULL
    colnames(y) <- y$samples$sample_id
    y[["samples"]][["group"]] <- y$samples$sample_type

    # remove some columns that came from the faciledataset
    axe.cols <- c("feature_type", "seqnames", "start", "end", "strand",
                  "effective_length", "source", "hdf5_index", "assay_type",
                  "assay")
    keep.cols <- !colnames(y[["genes"]]) %in% axe.cols
    y[["genes"]] <- y[["genes"]][, keep.cols, drop = FALSE]
  }

  fn.name <- paste0(".to_", class)
  fn <- try(getFunction(fn.name), silent = TRUE)
  if (is(fn, "try-error")) {
    stop("Conversion function not found: ", fn.name, call. = FALSE)
  }

  fn(y, ...)
}

#' @noRd
.to_DGEList <- function(x, ...) {
  x
}

#' @noRd
#' @importFrom stats model.matrix
.to_EList <- function(x, design = ~ group, ...) {
  limma::voom(x, model.matrix(design, data = x[["samples"]]))
}

#' @noRd
.to_ExpressionSet <- function(x, ...) {
  eset <- Biobase::ExpressionSet(x[["counts"]])
  eset <- Biobase::`pData<-`(eset, x[["samples"]])
  eset <- Biobase::`fData<-`(eset, x[["genes"]])
  eset
}

#' @noRd
.to_SummarizedExperiment <- function(x, ...) {
  SummarizedExperiment::SummarizedExperiment(
    list(counts = x[["counts"]]),
    rowData = x[["genes"]],
    colData = x[["samples"]])
}

#' @noRd
.to_DESeqDataSet <- function(x, design = ~ group, ...) {
  out <- DESeq2::DESeqDataSet(.to_SummarizedExperiment(x), design = design)
  DESeq2::estimateSizeFactors(out)
}

#' @noRd
.to_SingleCellExperiment <- function(x, ...) {
  stop("SingleCellExperiment support not yet implemented")
}

#' @noRd
.to_MultiAssayExperiment <- function(x, ...) {
  stop("MultiAssayExperiment support not yet implemented")
}
