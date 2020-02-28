#' @export
#' @noRd
biocbox.FacileBiocDataStore <- function(x, class = NULL,
                                        assay_name = default_assay(x),
                                        features = NULL, samples = NULL,
                                        custom_key = Sys.getenv("USER"), ...) {
  assert_choice(assay_name, assay_names(x))
  features.all <- features(x, assay_name = assay_name, ...)

  if (is.null(features)) {
    features <- features.all
  } else {
    if (is.data.frame(features)) {
      features <- features[["feature_id"]]
    }
    assert_character(features)
    features <- filter(features.all, feature_id %in% features)
  }

  x <- x[features[["feature_id"]],]

  if (!is.null(samples)) {
    samples. <- assert_sample_subset(samples)
    x <- x[, samples.[["sample_id"]]]
    x@facile[["assay_sample_info"]] <-
      lapply(x@facile[["assay_sample_info"]], semi_join, samples.,
             by = c("dataset", "sample_id"))
  }

  x
}
