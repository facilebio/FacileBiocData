if (!exists("Y")) Y <- NULL
if (!exists("dds")) dds <- example_bioc_data("DESeqDataSet", y = Y)
if (!exists("ddf")) ddf <- facilitate(dds)

test_that("facilitate.ExpressionSet works", {
  expect_s4_class(ddf, "FacileDESeqDataSet")
  expect_s4_class(ddf, "FacileBiocDataStore")
  expect_s4_class(ddf, "FacileDataStore")
  expect_s4_class(ddf, "DESeqDataSet")
  checkmate::expect_list(ifacile(ddf))

  covs <- samples(ddf) %>%
    with_sample_covariates("sex") %>%
    with_sample_covariates("sample_type")
  expect_s3_class(covs, "facile_frame")
  expect_equal(nrow(covs), unname(ncol(ddf)))
  checkmate::expect_factor(covs[["sex"]], c("m", "f"))
  checkmate::expect_factor(covs[["sample_type"]], c("normal", "tumor"))
})

test_that("counts(dds, normalized = TRUE) == fetch_assay_data(normalized = TRUE", {
  expected <- DESeq2::counts(dds, normalized = TRUE)
  # by default normalized = TRUE returns lognormal data. let's turn that off
  res <- fetch_assay_data(ddf, normalized = TRUE, log = FALSE, as.matrix = TRUE)
  expect_equal(rownames(res), rownames(expected))
  expect_equal(sub(".*__", "", colnames(res)), colnames(expected))
  expect_equal(res, expected, check.attributes = FALSE)
})

test_that("vst transformed data is retrievable from facilitated dds", {
  vsd <- SummarizedExperiment::assay(DESeq2::vst(dds, blind = TRUE))
  res <- fetch_assay_data(ddf, assay_name = "vst", as.matrix = TRUE)

  expect_equal(rownames(res), rownames(vsd))
  expect_equal(sub(".*__", "", colnames(res)), colnames(vsd))
  expect_equal(res, vsd, check.attributes = FALSE)
})

test_that("we know vst transformed data is lognorm & can be batch corrected", {
  res <- fetch_assay_data(ddf, assay_name = "vst", as.matrix = TRUE)
  resb <- fetch_assay_data(ddf, assay_name = "vst", batch = "sex", as.matrix = TRUE)
  diffs <- abs(res - resb)
  mean.diffs <- mean(as.vector(diffs))
  expect_gt(mean.diffs, 0.00)
  expect_lt(mean.diffs, 0.25)
})
