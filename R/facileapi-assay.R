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
assay_info.FacileBiocDataStore <- function(x, assay_name = NULL, ...) {
  anames <- assay_names(x)
  if (!is.null(assay_name)) {
    assert_choice(assay_name, assay_names)
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
.infer_assay_type <- function(x, amatrix, ...) {
  warning("TODO: .infer_assay_type needs serious improvement")
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

