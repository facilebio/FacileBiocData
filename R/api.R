#' Gets the FacileDataStore to play with S4 mojo(?)
#'
#' @noRd
setOldClass("FacileDataStore")


#' Root virtual class that signals object as a FacileDataStore
#'
#' All FacileBiocDataStore objects should inherit from this
#'
#' @noRd
#' @export
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

#' Retrieves names of assay elements in the container
#' @noRd
anames <- function(x, ...) {
  UseMethod("anames", x)
}
