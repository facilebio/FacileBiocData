context("ExpressionSet")

if (!exists("eset")) {
  bb <- loadNamespace("Biobase")
  eset <- local({
    efds <- FacileData::exampleFacileDataSet()
    y <- FacileData::as.DGEList(efds)

    y$samples$samid <- NULL
    colnames(y) <- y$samples$sample_id

    eset <- bb$ExpressionSet(y$counts)
    # pData(eset) <- y$samples
    # fData(eset) <- y$genes
    eset <- bb$`pData<-`(eset, y$samples)
    eset <- bb$`fData<-`(eset, y$genes)
    eset
  })
}

test_that("facilitate.ExpressionSet works", {
  f <- facilitate(eset)
  expect_s4_class(f, "FacileExpressionSet")
  expect_s4_class(f, "FacileBiocDataStore")
  expect_s4_class(f, "FacileDataStore")

  covs <- samples(f) %>%
    with_sample_covariates("sex") %>%
    with_sample_covariates("sample_type")
  expect_s3_class(covs, "facile_frame")
  expect_equal(nrow(covs), unname(ncol(f)))
  checkmate::expect_factor(covs[["sex"]], c("m", "f"))
  checkmate::expect_factor(covs[["sample_type"]], c("normal", "tumor"))
})
