#' @noRd
#' @export
samples.FacileBiocDataStore <- function(x, ...) {
  pdata(x) %>%
    as.tbl() %>%
    dplyr::select(dataset, sample_id) %>%
    as_facile_frame(x, .valid_sample_check = FALSE)
}
