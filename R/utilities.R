#' Preps the pdata(x) from a bioconductor container to be a legit facile-looking
#' sample frame
#'
#' This function is not exported on purpose.
#'
#' @param x a bioconductor container
#' @param sinfo the data.frame from `pdata(x)`
.init_pdata <- function(x, sinfo = pdata(x), ...) {
  stopifnot(is.data.frame(sinfo))
  if (is.null(sinfo[["sample_id"]])) sinfo[["sample_id"]] <- colnames(x)
  if (is.null(sinfo[["dataset"]])) sinfo[["dataset"]] <- "dataset"

  stopifnot(
    nrow(sinfo) == ncol(x),
    is.character(sinfo[["sample_id"]]),
    !any(duplicated(sinfo[["sample_id"]])),
    is.character(sinfo[["dataset"]]))

  rownames(sinfo) <- sinfo[["sample_id"]]
  sinfo
}

#' Preps the fdata() from bioconductor containers to look like a
#' facile_feature_frame
#' @noRd
#' @param x the bioconductor container
#' @param finfo the data.frame from `fdata(x)`
#' @param rename_fdata a character vector with `names(rename_fdata)` being
#'   the expected facile feature_info colnames()
#'   (.ie `colnames(.empte_feature_frame()`), the values are the colnames of
#'   `finfo()` that map to them
.init_fdata <- function(x, finfo = fdata(x), rename_fdata = NULL,
                        assay_name = NULL, assay_type = "unk_atype",
                        feature_type = "unk_ftype", ...) {
  stopifnot(is.data.frame(finfo))

  assay_names <- anames(x)
  if (is.null(assay_name)) {
    assay_name <- assay_names[1L]
  }
  assert_choice(assay_name, assay_names)

  assert_string(assay_type)
  assert_string(feature_type)

  ftmpl <- .empty_feature_frame()
  if (is.character(rename_fdata)) {
    valid <- names(rename_fdata) %in% colnames(ftmpl)
    rename_fdata <- rename_fdata[valid]
    if (length(rename_fdata)) {
      finfo <- rename(finfo, all_of(rename_fdata))
    }
  }

  # feature_id is really the only required thing, the rest we can make up
  if (is.character(finfo[["feature_id"]])) {
    if (any(duplicated(finfo[["feature_id"]]))) {
      warning("duplicated feature_id column detected, using rownames()")
      finfo[["feature_id"]] <- rownames(x)
    }
  } else {
    finfo[["feature_id"]] <- rownames(x)
  }
  if (any(duplicated(finfo[["feature_id"]]))) {
    stop("The feature 'feature_id' column is not unique")
  }
  if (!is.character(finfo[["name"]])) {
    # I wish I never renamed `name` to `symbol` in ye-old FacileData::as.DGEList
    if (is.character(finfo[["symbol"]])) {
      finfo[["name"]] <- finfo[["symbol"]]
    } else {
      finfo[["name"]] <- finfo[["feature_id"]]
    }
  }
  if (!is.character(finfo[["feature_type"]])) {
    finfo[["feature_type"]] <- feature_type
  }
  if (!is.character(finfo[["assay_type"]])) {
    finfo[["assay_type"]] <- assay_type
  }
  if (!is.character(finfo[["meta"]])) {
    finfo[["meta"]] <- "unk_meta"
  }
  for (fcol in colnames(ftmpl)) {
    expected_class <- class(ftmpl[[fcol]])[1L]
    if (!is(finfo[[fcol]], expected_class)) {
      stop("The '", fcol, "' feature column should be a subclass of type: ",
           expected_class)

    }
  }
  if (nrow(finfo) != nrow(x)) {
    stop("The number of rows in the feature data.frame is not the same as the ",
         "number of rows in the assay container")
  }
  finfo[["assay"]] <- assay_name
  finfo[["assay_type"]] <- assay_type
  finfo[["feature_type"]] <- feature_type
  rownames(finfo) <- finfo[["feature_id"]]
  finfo
}

.empty_feature_frame <- function(x = NULL) {
  out <- tibble(
    feature_id = character(),
    feature_type = character(),
    name = character(),
    meta = character())
  if (is(x, "FacileBiocDataStore")) {
    out <- as_facile_frame(out, x, .valid_sample_check = FALSE)
  }
  out
}
