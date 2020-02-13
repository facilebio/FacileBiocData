#' @include api.R
NULL

#' @export
#' @importFrom utils packageVersion
#' @importClassesFrom Biobase ExpressionSet
setClass("FacileExpressionSet",
         slots = c(facile = "list"),
         contains = c("FacileBiocDataStore", "ExpressionSet"),
         prototype = prototype(facile = list()))

#' @export
#' @noRd
#' @rdname facilitate
#' @method facilitate ExpressionSet
facilitate.ExpressionSet <- function(x, ...) {
  out <- new("FacileExpressionSet",
             experimentData = Biobase::experimentData(x),
             assayData = Biobase::assayData(x),
             phenoData = Biobase::phenoData(x),
             featureData = Biobase::featureData(x),
             anotation = Biobase::annotation(x),
             protocolData = Biobase::protocolData(x))
  out
}

# bioc data retrieval methods --------------------------------------------------

#' @noRd
fdata.ExpressionSet <- function(x, ...) {
  if (!requireNamespace("Biobase", quietly = TRUE)) {
    stop("Biobase package required, please install it.",
         call. = FALSE)
  }
  Biobase::fData(x)
}

#' @noRd
pdata.ExpressionSet <- function(x, ...) {
  if (!requireNamespace("Biobase", quietly = TRUE)) {
    stop("Biobase package required, please install it.",
         call. = FALSE)
  }
  Biobase::pData(x)
}

#' @noRd
adata.ExpressionSet <- function(x, name = "exprs", ...) {
  if (!requireNamespace("Biobase", quietly = TRUE)) {
    stop("Biobase package required, please install it.",
         call. = FALSE)
  }
  Biobase::assayDataElementy(x, name)
}
