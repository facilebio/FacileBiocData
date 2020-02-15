#' @noRd
#' @export
default_assay.FacileBiocDataStore <- function(x, ...) {
  out <- ifacile(x)[["default_assay"]]
  if (!test_string(out)) out <- assay_names(x, ...)[1L]
  out
}
