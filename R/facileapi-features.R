#' @noRd
#' @export
features.FacileBiocDataStore <- function(x, assay_name = default_assay(x),
                                         feature_type = NULL,
                                         feature_ids = NULL, ...,
                                         .developer = getOption("fbioc.developer", FALSE)) {
  assert_choice(assay_name, assay_names(x))
  if (is(x, "MultiAssayExperiment")) {
    stop("MUltiAssayExperiment support not yet implemented")
    # This is the only data type that has features from multiple universes,
    # otherwise the bioc container only has info for one type of feature.
    null.aname <- onull.aname <- is.null(assay_name)
    null.ftype <- onull.ftype <- is.null(feature_type)

    # Is the user asking for feature information from the features measured on
    # a given assay, or for all features of a given feature_type.
    if (null.aname && null.ftype) {
      assay_name <- default_assay(x)
      null.aname <- FALSE
    }
    if (!xor(null.aname, null.ftype)) {
      stop("Must specify feature information for EITHER assay_name or ",
           "feature_type")
    }
    if (null.aname) {
      query_type <- "feature_type"
      query_value <- assert_choice(feature_type, feature_types(x))
    } else {
      query_type <- "assay_name"
    }
  } else {
    if (!is.null(assay_name) || !is.null(feature_type)) {
      if (.developer) {
        warning("`assay_name` and `feature_type` are ignored for, '",
                class(x)[1L], "' data containers")
      }
    }
    out <- mutate(fdata(x), assay = assay_name)
  }

  if (!is.null(feature_ids)) {
    assert_character(feature_ids)
    out <- filter(out, .data$feature_id %in% feature_ids)
  }

  as_tibble(out)
}
