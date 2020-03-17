context("SummarizedExperiment")

if (!exists("Y")) Y <- NULL
if (!exists("dds")) dds <- example_bioc_data("DESeqDataSet", y = Y)

test_that("facilitate.ExpressionSet works", {
  f <- facilitate(dds)
  expect_s4_class(f, "FacileDESeqDataSet")
  expect_s4_class(f, "FacileBiocDataStore")
  expect_s4_class(f, "FacileDataStore")
  expect_s4_class(f, "DESeqDataSet")
  checkmate::expect_list(ifacile(f))

  covs <- samples(f) %>%
    with_sample_covariates("sex") %>%
    with_sample_covariates("sample_type")
  expect_s3_class(covs, "facile_frame")
  expect_equal(nrow(covs), unname(ncol(f)))
  checkmate::expect_factor(covs[["sex"]], c("m", "f"))
  checkmate::expect_factor(covs[["sample_type"]], c("normal", "tumor"))
})

# TODO: Add tests for fetching normalized assay data (vst, cpm)
