context("ExpressionSet")

if (!exists("Y")) Y <- NULL
if (!exists("eset")) eset <- example_bioc_data("ExpressionSet", y = Y)

test_that("facilitate.ExpressionSet works", {
  f <- facilitate(eset)
  expect_s4_class(f, "FacileExpressionSet")
  expect_s4_class(f, "FacileBiocDataStore")
  expect_s4_class(f, "FacileDataStore")
  expect_s4_class(f, "ExpressionSet")
  checkmate::expect_list(ifacile(f))

  covs <- samples(f) |>
    with_sample_covariates("sex") |>
    with_sample_covariates("sample_type")
  expect_s3_class(covs, "facile_frame")
  expect_equal(nrow(covs), unname(ncol(f)))
  checkmate::expect_factor(covs[["sex"]], c("m", "f"))
  checkmate::expect_factor(covs[["sample_type"]], c("normal", "tumor"))
})
