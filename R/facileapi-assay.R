#' @noRd
#' @export
assay_sample_info.FacileBiocDataStore <- function(
    x, assay_name, samples = NULL, ...,
    .developer = getOption("fbioc.developer", FALSE)) {
  assert_choice(assay_name, assay_names(x))

  if (is.null(samples)) {
    samples <- samples(x)
  } else {
    assert_sample_subset(samples, x)
    samples <- distinct(samples, .data[["dataset"]], .data[["sample_id"]],
                        .keep_all = TRUE)
  }

  samples <- collect(samples, n = Inf)
  assay.samples <- ifacile(x)[["assay_sample_info"]][[assay_name]]
  if (is.null(assay.samples)) {
    if (.developer) {
      warning("The '", assay_name, "' assay was added after call to ",
              "facilitate(), assay_sample_info data is missing")
    }
    assay.samples <- samples(x)
  }

  assay.samples |>
    semi_join(samples, by = c("dataset", "sample_id")) |>
    mutate(assay = .env$assay_name, .after = "dataset")
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
#' @export
#' @examples
#' yf <- example_bioc_data("DGEList") |> facilitate()
fetch_assay_data.FacileBiocDataStore <- function(
    x, features = NULL, samples = NULL, assay_name = default_assay(x),
    normalized = FALSE, batch = NULL, main = NULL, as.matrix = FALSE, ...,
    aggregate = FALSE, aggregate.by = "ewm", verbose = FALSE) {
  assert_flag(as.matrix)
  assert_flag(normalized)
  aggregate.by <- match.arg(tolower(aggregate.by), c("ewm", "zscore"))

  ainfo <- assay_info(x, assay_name)
  if (is.null(samples)) {
    samples <- collect(samples(x), n = Inf)
  } else {
    assert_sample_subset(samples)
    samples <- semi_join(samples, collect(samples(x), n = Inf),
                         by = c("dataset", "sample_id"))
  }

  if (is.null(assay_name)) assay_name <- default_assay(x)
  assert_choice(assay_name, assay_names(x))

  features.all <- features(x, assay_name)

  if (is.null(features)) {
    features <- features.all
  } else {
    if (is.data.frame(features)) {
      features <- features[["feature_id"]]
    }
    if (is.factor(features)) features <- as.character(features)
    if (is.character(features)) {
      features <- filter(features.all, .data$feature_id %in% .env$features)
    }
    stopifnot(is(features, 'tbl') || is(features, 'data.frame'))
    if (!'assay' %in% colnames(features) || !is.character(features$assay)) {
      features <- collect(features, n = Inf)
      features[["assay"]] <- assay_name
    }
    assert_assay_feature_descriptor(features)
  }

  features <- distinct(features, .data$feature_id, .keep_all = TRUE)

  samples[["samid"]] <- paste(samples$dataset, samples$sample_id, sep = "__")
  adat.all <- adata(x, assay_name)[, samples[["samid"]], drop = FALSE]
  adat <- adat.all[features[["feature_id"]],, drop = FALSE]
  features[["assay_type"]] <- ainfo[["assay_type"]]

  if (!is.null(batch)) {
    if (!isTRUE(normalized)) {
      warning("`batch` parameter specified, setting `normalized` to `TRUE`")
    }
    normalized <- TRUE
  }

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
    adat <- normalize_assay_data(adat, features, samples, batch = batch,
                                 main = main, verbose = verbose, ...)
  }

  pdat <- pdata(x)
  samples[["samid"]] <- NULL

  aggregated <- NULL

  if (isTRUE(aggregate)) {
    if (aggregate.by == "ewm") {
      aggregated <- sparrow::eigenWeightedMean(adat, ...)
    } else if (aggregate.by == "zscore") {
      aggregated <- sparrow::zScore(adat, ...)
    } else {
      stop("Unknown aggregation method: ", aggregate.by)
    }
    adat <- matrix(
      aggregated$score,
      nrow = 1L,
      dimnames = list("score", names(aggregated$score)))
  }

  if (!as.matrix) {
    atype <- ainfo[["assay_type"]]
    ftype <- ainfo[["feature_type"]]
    vals <- FacileData:::.melt.assay.matrix(adat, assay_name, atype, ftype,
                                            features)
    vals <- as_tibble(vals)
    if (isTRUE(aggregate)) {
      vals <- mutate(vals,
                     feature_type = "aggregated",
                     feature_id = "aggregated",
                     feature_name = "aggregated")
    }
    vals <- samples |>
      left_join(
        vals,
        by = c("dataset", "sample_id"),
        suffix = c("", ".dropme..")) |>
      select(-ends_with(".dropme.."))
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

  anames <- assay_names(x, ...)
  if (!is.null(assay_name)) {
    assert_choice(assay_name, anames)
    anames <- assay_name
  }

  assay_info. <- ifacile(x)[["assay_info"]]
  # adding names(assay_info.) because we add "ghost" assays for some type
  # of containers, ie. the DESeqDataSet has a "normcounts" ghost assay
  # anames <- unique(anames, names(assay_info.))
  # maybe not

  ainfo <- lapply(anames, function(aname) {
    adat <- adata(x, aname)
    cached <- assay_info.[[aname]]

    ai <- list(
      assay = aname,
      assay_type = cached[["assay_type"]],
      feature_type = cached[["feature_type"]],
      description = cached[["description"]],
      nfeatures = nrow(adat),
      storage_mode = class(adat[1L])[1L])

    if (is.null(ai[["feature_type"]])) {
      # finfo <- suppressWarnings(FacileData::infer_feature_type(rownames(adat)))
      finfo <- fdata(x)
      ftype <- finfo[["id_type"]]
      if (is.null(ai[["feature_type"]])) {
        if (length(unique(ftype)) == 1L) {
          ftype <- ftype[1L]
        } else {
          # warning("Mixed feature_types in assay: ", aname, immediate. = TRUE)
          ftype <- "mixed"
        }
        ai[["feature_type"]] <- ftype
      }
    }

    if (is.null(ai[["assay_type"]])) {
      ai[["assay_type"]] <- .infer_assay_type(x, assay_name = aname, ...)
    }

    if (is.null(ai[["description"]])) {
      ai[["description"]] <- paste("assay data from '", aname, "'", sep = "")
    }

    as_tibble(ai)
  })

  bind_rows(ainfo)
}

# Internal Helpers -------------------------------------------------------------

#' @noRd
#' @param x the BiocDataContainer the assay came from
#' @param amatrix an assay matrix
.infer_assay_type <- function(x, assay_name, assay_type = NULL, ...,
                              .developer = getOption("fbioc.developer", FALSE)) {
  if (test_string(assay_type)) return(assay_type)

  if (.developer) {
    warning("TODO: .infer_assay_type needs serious improvement")
  }
  amatrix <- assert_matrix(adata(x, assay_name))
  atype <- NULL

  if (test_multi_class(x, .bioc_assume_count_container) &&
      assay_name == default_assay(x)) {
    # NOTE: this should be "counts"?
    return("rnaseq")
  }

  asummary <- summary(as.vector(amatrix))
  if (asummary["Min."] < 0 & asummary["Max."] < 20) return("lognorm")

  return("raw")
}

