#' Gets the FacileDataStore to play with S4 mojo(?)
#' @noRd
#' @export
setOldClass("FacileDataStore")


#' Root virtual class that signals object as a FacileDataStore
#'
#' All FacileBiocDataStore objects should inherit from this
#'
#' @noRd
#' @export
setClass("FacileBiocDataStore",
         # representation("VIRTUAL"),
         contains = c("FacileDataStore", "VIRTUAL"))

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
