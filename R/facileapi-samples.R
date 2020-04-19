#' @noRd
#' @export
samples.FacileBiocDataStore <- function(x, ...) {
  pdata(x) %>%
    as.tbl() %>%
    dplyr::select(.data$dataset, .data$sample_id) %>%
    as_facile_frame(x, .valid_sample_check = FALSE)
}

#' @noRd
#' @export
filter_samples.FacileBiocDataStore <- function(x, ..., samples. = samples(x),
                                               custom_key = Sys.getenv("USER"),
                                               with_covariates = FALSE) {
  force(samples.)
  assert_sample_subset(samples.)

  out <- pdata(x) %>%
    semi_join(samples., by = c("dataset", "sample_id")) %>%
    filter(...)
  if (!with_covariates) {
    out <- select(out, .data$dataset, .data$sample_id)
  }
  if (nrow(out) == 0L) {
    warning("All samples have been filtered out", immediate. = TRUE)
  }
  as_facile_frame(out, x, .valid_sample_check = FALSE)
}
