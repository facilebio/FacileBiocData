#' @export
#' @noRd
biocbox.FacileBiocDataStore <- function(x, class = NULL,
                                        assay_name = default_assay(x),
                                        features = NULL, samples = NULL,
                                        custom_key = Sys.getenv("USER"), ...) {
  assert_choice(assay_name, assay_names(x))
  if (!is.null(class)) {
    if (class == "list") {
      out <- biocbox(samples(x), class = "list", assay_name = assay_name,
                     features = features, custom_key = custom_key, ...)
      return(out)
    } else {
      if (!is(x, class)) {
        stop("Conversion from one Bioc data container to another is not, ",
             "supported. Leave class parameter as `NULL`")
      }
    }
  }

  features.all <- features(x, assay_name = assay_name, ...)

  if (!is.null(features)) {
    if (is.data.frame(features)) features <- features[["feature_id"]]
    assert_character(features)
    if (!test_subset(features, features.all[["feature_id"]])) {
      stop("features requested that do not exist in the FacileBiocDataStore")
    }
    # FIXME: w ith MultiAssayExperiment, subsetting features like this won't work
    x <- x[features,]
  }

  if (!is.null(samples)) {
    samples. <- assert_sample_subset(samples, fds = x)
    if (nrow(samples.) == 0) {
      stop("zero samples selected, can not return an empty data container")
    }
    # retrieve indices of requested samples, let's return the dataset in the
    # same order that samples came in
    idx <- match(
      with(samples., paste0(dataset, sample_id)),
      with(samples(x), paste0(dataset, sample_id)))

    x <- x[, idx]
    x@facile[["assay_sample_info"]] <-
      lapply(x@facile[["assay_sample_info"]], semi_join, samples.,
             by = c("dataset", "sample_id"))
  }

  x
}
