#' Gets the FacileDataStore to play with S4 mojo(?)
#'
#' @noRd
setOldClass("FacileDataStore")

#' Bioconductor classes that we assume/assert provide count data, like
#' rnaseq.
#' @noRd
.bioc_assume_count_container <- c(
  "DGEList", "DESeqDataSet", "SingleCellExperiment")

#' Immerse bioconductor assay containers into the facile.bio ecosystem.
#'
#' Bioconductor assay containers, like a DGEList, DESeqDataSet,
#' SummarizedExperiment, etc. can be used within the facie.bio ecosystem by
#' invoking the `facilitate()` method on them. This will return a `Facile*`
#' subclass of the container itself.
#'
#' For instance, `facilitate(DGEList)` will return a `FacileDGEList`, which can
#' be used a "normal" DGEList in all the same ways, but is also wrapped with
#' the facile api api and can be used by methods withing the `FacileAnalysis`,
#' for instance.
#'
#' These classes are also all subclass of the abstract `FacileBiocDataStore`
#' virtual class.
#'
#' @export
#' @aliases facilitate
#' @aliases FacileBiocDataStore
#'
#' @param assay_type A string that indicates the type of assay stored in the
#'   primary assay of the container. For some assay containers, like
#'   `DESeqDataSet`, `DGEList`, and `SingleCellExperiment`, we can assume the
#'   default value for this to be `"rnaseq"`. For the rest, we assume it's
#'   `"lognorm"`.
#' @param feature_type A string that indicates the type of features identifiers
#'   the assay containers is using. Default is `"infer"` to try to guess, but
#'   this is not the most accurate.
#' @param organism the organism the dataset is for (Homo sapiens, Mus musculus,
#'   etc.)
#' @param verbose make some noise
setClass("FacileBiocDataStore",
         contains = c("FacileDataStore", "VIRTUAL"),
         slots = c(facile = "list"),
         prototype = prototype(facile = list()))

#' Retrieve the internal facile list data structure from a data container
#'
#' @noRd
ifacile <- function(x, ...) {
  UseMethod("ifacile", x)
}

#'@noRd
ifacile.FacileBiocDataStore <- function(x, ...) {
  x@facile
}

#
# Data retrieval methods for bioconductor objects ------------------------------
#
# These functions are intentionally not exported

#' Retrieves feature-level metadata from a bioc container
#' @noRd
fdata <- function(x, ...) {
  UseMethod("fdata", x)
}

#' Retrieves sample-level metadata from a bioc container
#' @noRd
pdata <- function(x, ...) {
  UseMethod("pdata", x)
}

#' Retrieves assay-level data from a bioc container
#' @noRd
adata <- function(x, assay = NULL, ...) {
  UseMethod("adata", x)
}

#' some bioc containers allow assay matrices w/ NULL rownames(), like the
#' data in the airway SE. Currently we assume the rownames of the assay
#' matrix are feature_ids, so making sure that's the case.
#' @noRd
.cleanup_adata <- function(x, adat, ...) {
  ids <- rownames(adat)
  # the rownames of the airway count matrix are NULL, so we have to be
  # careful
  if (is.null(ids)) {
    if (nrow(adat) != nrow(x)) {
      stop("Feature space of assay matrix != feature space of data ",
           "container, is this a MultiAssayExperiment?\n  ",
           "https://github.com/facilebio/FacileBiocData/issues/1")
    }
    rownames(adat) <- rownames(x)
  }

  bad.rownames <- is.null(rownames(adat))
  if (bad.rownames) {
    stop("rownames of assay matrix could not be set to feature ids")
  }

  fids.all <- fdata(x, ...)[["feature_id"]]
  if (is.character(fids.all) && length(setdiff(rownames(adat), fids.all))) {
    stop("rownames of assay matrix not found in feature_id universe")
  }

  adat
}

