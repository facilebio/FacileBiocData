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
