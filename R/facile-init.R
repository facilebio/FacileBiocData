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

  rownames(sinfo) <- paste(sinfo[["dataset"]], sinfo[["sample_id"]], sep = "__")
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
                        feature_type = "infer", feature_source = "unknown",
                        ...) {
  stopifnot(is.data.frame(finfo))

  # assay_names <- assay_names(x)
  # if (is.null(assay_name)) {
  #   assay_name <- assay_names[1L]
  # }
  # assert_choice(assay_name, assay_names)
  #
  # assert_string(assay_type)
  # assert_string(feature_type)

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

  if (!is.character(finfo[["symbol"]])) {
    finfo[["symbol"]] <- finfo[["name"]]
  }

  if (!is.character(finfo[["feature_type"]]) ||
      !(missing(feature_type) && feature_type[1L] == "infer")) {
    ftype <- FacileData::infer_feature_type(finfo[["feature_id"]])
    stopifnot(all.equal(ftype[["id"]], finfo[["feature_id"]]))
    ftype <- ftype[["id_type"]]
    if (length(unique(ftype)) > 1L) {
      warning("Mixed feature_types in assay fdata()", immediate. = TRUE)
    }
    finfo[["feature_type"]] <- ftype
  }

  if (!is.character(finfo[["meta"]])) {
    finfo[["meta"]] <- "unk_meta"
  }

  if (!is.character(finfo[["source"]])) {
    finfo[["source"]] <- feature_source
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
  rownames(finfo) <- finfo[["feature_id"]]
  finfo
}

#' Initialize some sample-level assay information, like normfactor/sizefactor
#' etc.
#'
#' This was first developed to get the libsize and normfactor columns onto a
#' samples frame for fetch_assay_data(..., normalized = TRUE) to work.
#'
#' **NOTE**: currently any rnaseq assay will return normalized value using
#' edgeR::cpm mojo, even a DESeqDataSet
#'
#' Each assay will be given its own dataset,sample table where assay-level
#' metadata can be stored.
#' @noRd
#' @importFrom Matrix colSums
#' @importFrom edgeR calcNormFactors
.init_assay_sample_info <- function(x, ...) {
  assert_class(x, "FacileBiocDataStore")
  # this is only to support "ghost" assays like "cpm" in DESeqDataSet
  ainfo <- assay_info(x)
  pdat <- pdata(x)
  samples. <- samples(x)
  out <- lapply(seq(nrow(ainfo)), function(i) {
    info <- ainfo[i,,drop = FALSE]
    adat <- adata(x, info[["assay"]])

    if (is(x, "DGEList")) {
      samples. <- mutate(samples., libsize = pdat[["lib.size"]],
                         normfactor = pdat[["norm.factors"]])
      assert_numeric(samples.[["libsize"]])
      assert_numeric(samples.[["normfactor"]])
    }
    if (info[["assay_type"]] == "rnaseq") {
      lsize <- samples.[["libsize"]]
      nfctr <- samples.[["normfactor"]]
      if (!is.numeric(lsize)) {
        samples.[["libsize"]] <- Matrix::colSums(adat)
      }
      if (!is.numeric(nfctr) || all(nfctr == 1)) {
        samples.[["normfactor"]] <- calcNormFactors(adat, samples.[["libsize"]])
      }
    }
    samples.
  })
  names(out) <- ainfo[["assay"]]
  out
}

#' @noRd
.empty_feature_frame <- function(x = NULL) {
  out <- tibble(
    feature_id = character(),
    feature_type = character(),
    name = character(),
    meta = character(),
    source = character())
  if (is(x, "FacileBiocDataStore")) {
    out <- as_facile_frame(out, x, .valid_sample_check = FALSE)
  }
  out
}
