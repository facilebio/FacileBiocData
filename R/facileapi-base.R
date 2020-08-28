#' @noRd
#' @export
name.FacileBiocDataStore <- function(x, ...) {
  name <- ifacile(x)[["name"]]
  if (is.null(name)) {
    name <- paste(colnames(x), collapse = ":")
  }
  name
}

#' @noRd
#' @export
organism.FacileBiocDataStore <- function(x, ...) {
  out <- ifacile(x)[["organism"]]
  if (is.null(out)) out <- "unknown"
  out
}
