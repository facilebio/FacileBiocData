#' @noRd
#' @export
organism.FacileBiocDataStore <- function(x, ...) {
  out <- ifacile(x)[["organism"]]
  if (is.null(out)) out <- "unknown"
  out
}
