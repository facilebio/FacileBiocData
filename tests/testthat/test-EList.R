context("EList")

if (!exists("Y")) Y <- NULL
if (!exists("elist")) elist <- example_bioc_data("EList", y = Y)

test_that("facilitate.EList works", {
  f <- facilitate(elist)

  expect_setequal(names(f), names(elist))
  expect_s4_class(f, "FacileEList")
  expect_s4_class(f, "FacileBiocDataStore")
  expect_s4_class(f, "FacileDataStore")
  expect_s4_class(f, "EList")
  checkmate::expect_list(ifacile(f))

  covs <- samples(f) |>
    with_sample_covariates("sex") |>
    with_sample_covariates("sample_type")
  expect_s3_class(covs, "facile_frame")
  expect_equal(nrow(covs), unname(ncol(f)))
  checkmate::expect_factor(covs[["sex"]], c("m", "f"))
  checkmate::expect_factor(covs[["sample_type"]], c("normal", "tumor"))

  expect_s4_class(fds(covs), "FacileEList")
})
