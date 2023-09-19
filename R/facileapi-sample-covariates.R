#' This function isn't optimized much. It always re-creates an EAVtable from
#' the pdata(x) data.frame. The only optimization is that it first reduces
#' the rows and columns of `pdata(x)` to match the `samples` (rows) and
#' `covariates` (columns) that were called for.
#'
#' The respective EAV covariate definition is always reconstructed in
#' `covariate_definitions.FacileBiocDataStore`
#'
#' By only defininng `fetch_sample_covariates.FacileBiocDataStore`, we kill
#' don't also have to define `with_sample_covariates.FacileBiocDataStore`,
#' however when `with_` first calls this, we have to melt and recast again.
#'
#' @noRd
#' @export
fetch_sample_covariates.FacileBiocDataStore <- function(
  x, covariates = NULL, samples = NULL,
  custom_key = Sys.getenv("USER"), with_source = FALSE, ...) {

  sinfo <- pdata(x)
  if (is.null(samples)) {
    samples <- samples(x)
  } else {
    assert_sample_subset(samples)
    sinfo <- semi_join(sinfo, samples, by = c("dataset", "sample_id"))
  }

  all.covs <- setdiff(colnames(sinfo), c("dataset", "sample_id"))
  if (is.null(covariates)) {
    covariates <- all.covs
  } else {
    assert_character(covariates)
    covariates <- setdiff(covariates, c("dataset", "sample_id"))
    bad.covs <- setdiff(covariates, all.covs)
    if (length(bad.covs)) {
      warning("Unknown sample covariates: ", paste(bad.covs, collapse = ","))
      covariates <- intersect(covariates, all.covs)
    }
    sinfo <- select(sinfo, "dataset", "sample_id", {{covariates}})
  }

  if (nrow(samples) == 0 || length(all.covs) == 0L) {
    return(.empty_sample_covariates(x))
  }

  out <- as.EAVtable(sinfo)
  as_facile_frame(out, x, .valid_sample_check = FALSE)
}

#' @noRd
#' @export
covariate_definitions.FacileBiocDataStore <- function(x, as.list = TRUE, ...) {
  # eav <- as.EAVtable(pdata(x))
  # out <- attr(eav, "covariate_def")
  # out <- ifacile(x)[["covariate_def"]]
  out <- eav_metadata_create(pdata(x))
  if (!as.list) {
    out <- lapply(names(out), function(name) {
      i <- out[[name]]
      lvls <- unique(i$levels)
      is.factor <- !is.null(lvls)
      lbl <- if (is.null(i$label)) name else i$label
      tibble(variable=name, type=i$type, class=i$class, label=i$label,
             is_factor=is.factor, levels=list(lvls), description=i$description)
    })
    out <- bind_rows(out)
  }
  class(out) <- c('CovariateDefinitions', class(out))
  set_fds(out, x)
}

#' @noRd
#' @export
fetch_custom_sample_covariates.FacileBiocDataStore <- function(
  x, covariates = NULL, samples = NULL, custom_key = Sys.getenv("USER"),
  file.prefix = "facile", ...) {
  .empty_sample_covariates(x)
}

#' @noRd
.empty_sample_covariates <- function(x = NULL) {
  out <- tibble(
    dataset = character(),
    sample_id = character(),
    variable = character(),
    value = character(),
    class = character(),
    type = character(),
    date_entered = integer())
  if (is(x, "FacileBiocDataStore")) {
    out <- as_facile_frame(out, x, .valid_sample_check = FALSE)
  }
  out
}
