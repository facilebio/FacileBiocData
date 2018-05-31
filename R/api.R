#' Creates FacileData sublcasses of Bioconductor assay containers
#'
#' This function creates a subclass of given Bioconductor assay container that
#' implements the FacileData API. This is in contrast with the
#' [FacileData::as.FacileDataSet()] function, which literally converts an input
#' object into a completely new `FacileDataSet`.
#'
#' @export
#' @param x A bioconductor assay container
#' @param ... parameters to tweak container-specific funcitonality
#' @return A facile-subclass of `x` that implements teh FacileData API
facilitate <- function(x, ...) {
  UseMethod("facilitate", x)
}

