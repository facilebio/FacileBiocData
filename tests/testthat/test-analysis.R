context("FacileAnalysis")

# suppressWarnings(library(FacileAnalysis))
devtools::load_all(".")
devtools::load_all("/Users/lianoglou/workspace/facilebio/public/packages/FacileAnalysis")

if (!exists("dlist")) {
  dlist <- example_bioc_data("DGEList")
}

test_that("flm_def defines t-test and anova models on FacileBiocDataStore", {
  f <- facilitate(dlist)
  des.ttest <- flm_def(f, "sample_type", "tumor", "normal", batch = "sex")
  expect_s3_class(des.ttest, "FacileTtestModelDefinition")

  des.anova <- flm_def(f, covariate = "stage", batch = "sex")
  expect_s3_class(des.anova, "FacileAnovaModelDefinition")
})

test_that("biocbox can be constructed from a FacileLinearModelDefinition", {
  f <- facilitate(dlist)
  des.ttest <- flm_def(f, covariate = "sample_type", "tumor", "normal",
                       batch = "sex")
  des <- design(des.ttest)
  expect_equal(
    rownames(des),
    with(samples(f), paste(dataset, sample_id, sep = "__")))

  bb <- biocbox(des.ttest)
  expect_s3_class(des.ttest, "FacileTtestModelDefinition")
  des.anova <- flm_def(f, covariate = "stage", batch = "sex")
  expect_s3_class(des.anova, "FacileAnovaModelDefinition")
})


test_that("flm_def defines models on facile_from from a FacileBiocDataStore", {
  f <- facilitate(dlist)
  blca <- filter_samples(f, indication == "BLCA")

  des.ttest <- flm_def(blca, covariate = "sample_type", "tumor", "normal",
                       batch = "sex")
  expect_s3_class(des.ttest, "FacileTtestModelDefinition")

  des.anova <- flm_def(blca, covariate = "stage", batch = "sex")
  expect_s3_class(des.anova, "FacileAnovaModelDefinition")

  # The `blca` sample didn't have any covariates, but thy get spanked on within
  # flm_def so we can keep track of what the covariates were when tested.
  # Retrieving those covariates back to blca should make it equiavlent again
  tsamples <- with_sample_covariates(blca, c("sample_type", "sex"))
  asamples <- with_sample_covariates(blca, c("sex", "stage"))
  expect_equal(samples(des.ttest), tsamples)
  expect_equal(samples(des.anova), asamples)
})

test_that("fdge can voom", {
  f <- facilitate(dlist)
  blca <- filter_samples(f, indication == "BLCA")
  des.ttest <- flm_def(blca, covariate = "sample_type", "tumor", "normal",
                       batch = "sex")
  vmf <- fdge(des.ttest, method = "voom")
})
