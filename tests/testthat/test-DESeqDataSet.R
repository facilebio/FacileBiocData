context("SummarizedExperiment")

if (!exists("dds")) {
  ns <- loadNamespace("DESeq2")
  dds <- local({
    efds <- FacileData::exampleFacileDataSet()
    y <- FacileData::as.DGEList(efds)

    y$samples$samid <- NULL
    colnames(y) <- y$samples$sample_id

    se <- SummarizedExperiment::SummarizedExperiment(
      y$counts, rowData = y$genes, colData = y$samples)
    ns$DESeqDataSet(se, design = ~ sample_type)
  })
}

test_that("facilitate.ExpressionSet works", {
  f <- facilitate(dds)
  expect_s4_class(f, "FacileDESeqDataSet")
  expect_s4_class(f, "FacileBiocDataStore")
  expect_s4_class(f, "FacileDataStore")
  expect_s4_class(f, "DESeqDataSet")

  covs <- samples(f) %>%
    with_sample_covariates("sex") %>%
    with_sample_covariates("sample_type")
  expect_s3_class(covs, "facile_frame")
  expect_equal(nrow(covs), unname(ncol(f)))
  checkmate::expect_factor(covs[["sex"]], c("m", "f"))
  checkmate::expect_factor(covs[["sample_type"]], c("normal", "tumor"))
})
