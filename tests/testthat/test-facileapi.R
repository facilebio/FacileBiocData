context("Facile API over facilitated containers")


.classes <- c("DESeqDataSet", "DGEList", "EList", "ExpressionSet",
              "SummarizedExperiment")
if (!exists("Y")) Y <- example_bioc_data("DGEList")

BIOC <- sapply(.classes, example_bioc_data, Y = Y, simplify = FALSE)

test_that("assay_names is defined for base Bioconductor classes", {
  for (bclass in names(BIOC)) {
    obj <- BIOC[[bclass]]
    fn.name <- paste0("assay_names.", class(obj)[1L])
    checkmate::expect_function(getFunction(fn.name), info = bclass)
    anames <- assay_names(obj)
    checkmate::expect_character(anames, min.len = 1L, info = bclass)
  }
})

test_that("assay_names returns names of assays in the facilitated container", {
  for (bclass in names(BIOC)) {
    obj <- BIOC[[bclass]]
    f <- facilitate(obj)
    checkmate::expect_set_equal(assay_names(f), assay_names(obj), info = bclass)
  }
})

test_that("assay_info returns legit metadata for all containers", {
  expected.cols <- c(
    assay = "character",
    assay_type = "character",
    feature_type = "character",
    description = "character",
    nfeatures = "integer",
    storage_mode = "character")

  for (bclass in names(BIOC)) {
    obj <- BIOC[[bclass]]
    f <- facilitate(obj)
    ainfo <- expect_warning(assay_info(f), "improve")
    expect_s3_class(ainfo, "data.frame")

    # check each column is of expected data type
    for (cname in names(expected.cols)) {
      ctype <- expected.cols[cname]
      info <- sprintf("column '%s' is not type '%s' from container '%s'",
                      cname, ctype, bclass)
      expect_is(ainfo[[cname]], ctype, info = info)
    }

    expect_equal(ainfo[["feature_type"]][1L], "entrez", info = blcass)
    expect_equal(ainfo[["nfeatures"]][1L], unname(nrow(obj)), info = bclass)
    checkmate::expect_set_equal(ainfo[["assay"]], assay_names(f), info = bclass)
  }
})
