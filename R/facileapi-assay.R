#' @noRd
#' @export
assay_sample_info.FacileBiocDataStore <- function(x, assay_name, samples = NULL,
                                                  ...) {
  assert_choice(assay_name, assay_names(x))
  if (is.null(samples)) {
    samples <- samples(x)
  } else {
    assert_sample_subset(samples, x)
    samples <- distinct(samples, dataset, sample_id, .keep_all = TRUE)
  }

  samples <- collect(samples, n = Inf)
  semi_join(ifacile(x)[["assay_sample_info"]][[assay_name]], samples,
            by = c("dataset", "sample_id"))
}

#' Assay data retrieval from a FacileBiocDataStore
#'
#' required output columns:
#' * dataset [chr]
#' * sample_id [chr]
#' * assay [chr]
#' * assay_type [chr]
#' * feature_type [chr]
#' * feature_id [chr]
#' * feature_name [chr]
#' * value [int]
#'
#' @noRd
#' @examples
#' yf <- example_bioc_data("DGEList") %>% facilitate()
fetch_assay_data.FacileBiocDataStore <- function(
    x, features, samples = NULL, assay_name = default_assay(x),
    normalized = FALSE, batch = NULL, main = NULL, as.matrix = FALSE,
    ..., aggregate = FALSE, aggregate.by= "ewm", verbose = FALSE) {
  assert_flag(as.matrix)
  assert_flag(normalized)
  if (!normalized && !is.null(batch)) {
    warning("No batch correction will happen when normalized = FALSE")
  }
  ainfo <- assay_info(x, assay_name)
  if (is.null(samples)) {
    samples <- collect(samples(x), n = Inf)
  } else {
    assert_sample_subset(samples)
    samples <- semi_join(samples, collect(samples(x), n = Inf),
                         by = c("dataset", "sample_id"))
  }

  if (!is.null(assay_name) || is.character(assay_name)) {
    assert_string(assay_name)
    assert_choice(assay_name, assay_names(x))
  }

  if (missing(features) || is.null(features)) {
    assert_string(assay_name)
    features <- FacileBiocData::features(x, assay_name)
  } else {
    if (is.character(features)) {
      features <- tibble(feature_id=features, assay=assay_name)
    }
    stopifnot(is(features, 'tbl') || is(features, 'data.frame'))
    if (!'assay' %in% colnames(features) || !is.character(features$assay)) {
      features <- collect(features, n = Inf)
      features[["assay"]] <- assay_name
    }
    assert_assay_feature_descriptor(features)
  }

  features <- distinct(features, feature_id, .keep_all = TRUE)

  adat.all <- adata(x, assay_name)[, samples[["sample_id"]], drop = FALSE]
  adat <- adat.all[features[["feature_id"]],, drop = FALSE]
  colnames(adat) <- with(samples, paste(dataset, sample_id, sep = "__"))
  features[["assay_type"]] <- ainfo[["assay_type"]]

  samples[["samid"]] <- colnames(adat)

  if (normalized) {
    # Adds sample-level assay data appropriate for whatever the assay is
    asinfo <- assay_sample_info(x, assay_name, samples)
    # If `samples` were passed in with any assay-level covariates, let those
    # override what is in the database
    custom.cols <- intersect(colnames(samples), colnames(asinfo))
    custom.cols <- setdiff(custom.cols, c("dataset", "sample_id"))
    if (length(custom.cols)) {
      asinfo <- asinfo[, !colnames(asinfo) %in% custom.cols]
    }
    if (length(setdiff(colnames(asinfo), c("dataset", "sample_id")))) {
      samples <- left_join(samples, asinfo, by = c("dataset", "sample_id"))
    }
    adat <- normalize_assay_data(adat, features, samples,
                                 batch = batch, main = main,
                                 verbose = verbose, ...)
  }

  pdat <- pdata(x)
  samples[["samid"]] <- NULL

  if (!as.matrix) {
    atype <- ainfo[["assay_type"]]
    ftype <- ainfo[["feature_type"]]
    vals <- FacileData:::.melt.assay.matrix(adat, assay_name, atype, ftype,
                                            features)
    vals <- as.tbl(vals)
    if (isTRUE(aggregate)) {
      vals <- mutate(vals,
                     feature_type = "aggregated",
                     feature_id = "aggregated",
                     feature_name = "aggregated")
    }
    vals <- left_join(samples, vals, by = c("dataset", "sample_id"))
  } else {
    vals <- adat
  }

  set_fds(vals, x)
}

#' @noRd
#' @export
default_assay.FacileBiocDataStore <- function(x, ...) {
  out <- ifacile(x)[["default_assay"]]
  if (!test_string(out)) out <- assay_names(x, ...)[1L]
  out
}

#' Required for FacileAnalysis::fdge
#'
#' @noRd
#' @export
assay_info.FacileBiocDataStore <- function(x, assay_name = NULL, ...,
                                           .developer = getOption("fbioc.developer", FALSE)) {
  if (.developer) {
    warning("TODO: assay_info.FacileBiocDataStore needs improvement")
  }

  anames <- assay_names(x)
  if (!is.null(assay_name)) {
    assert_choice(assay_name, anames)
    anames <- assay_name
  }

  ainfo <- lapply(anames, function(aname) {
    adat <- adata(x, aname)
    finfo <- FacileData::infer_feature_type(rownames(adat))
    ftype <- finfo[["id_type"]]

    if (length(unique(ftype)) == 1L) {
      ftype <- ftype[1L]
    } else {
      warning("Mixed feature_types in assay: ", aname, immediate. = TRUE)
      ftype <- "mixed"
    }

    tibble(
      assay = aname,
      assay_type = .infer_assay_type(x, adat, ...),
      feature_type = ftype,
      description = paste("assay data from '", aname, "'", sep = ""),
      nfeatures = nrow(adat),
      storage_mode = class(adat[1L])[1L])
  })

  bind_rows(ainfo)
}

# Internal Helpers -------------------------------------------------------------

#' @noRd
#' @param x a FacileBiocDataStore
.init_assay_info <- function(x, assay_type = "infer", feature_type = "infer",
                             ...) {

}

#' @noRd
#' @param x the BiocDataContainer the assay came from
#' @param amatrix an assay matrix
.infer_assay_type <- function(x, amatrix, ...,
                              .developer = getOption("fbioc.developer", FALSE)) {
  if (.developer) {
    warning("TODO: .infer_assay_type needs serious improvement")
  }
  assert_matrix(amatrix)
  atype <- NULL

  rnaseq.class <- c("DGEList", "DESeqDataSet", "SingleCellExperiment")
  if (test_multi_class(x, rnaseq.class)) {
    return("rnaseq")
  }

  asummary <- summary(as.vector(amatrix))
  if (asummary["Min."] < 0 & asummary["Max."] < 20) {
    return("lognorm")
  }
  if (asummary["Min."] >= 0) {
    return("rnaseq")
  }

  stop(".infer_assay_type This needs to be improved", call. = FALSE)
}

