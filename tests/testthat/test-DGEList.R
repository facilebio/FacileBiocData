context("DGEList")

library(edgeR)

if (!exists("dgelist")) {
  dlist <- local({
    efds <- FacileData::exampleFacileDataSet()
    y <- FacileData::as.DGEList(efds)
    y$samples$samid <- NULL
    colnames(y) <- y$samples$sample_id
    y
  })
}

test_that("facilitate.DGEList works", {
  f <- facilitate(dlist)
  expect_s4_class(f, "FacileDGEList")
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
