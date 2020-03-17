context("FacileAnalysis: Differential Expression")

.classes <- c("DESeqDataSet", "DGEList", "EList", "ExpressionSet",
              "SummarizedExperiment")
.rnaseq.class <- setdiff(.classes, "EList")

if (!exists("FDS")) FDS <- FacileData::exampleFacileDataSet()
if (!exists("Y")) Y <- example_bioc_data("DGEList", efds = FDS)

BIOC <- sapply(.classes, example_bioc_data, Y = Y, simplify = FALSE)
FBIOC <- lapply(BIOC, function(b) {
  assay_type <- if (is(b, "EList")) "lognorm" else "rnaseq"
  facilitate(b, assay_type = assay_type, run_vst = FALSE)
})

suppressPackageStartupMessages(suppressWarnings(library(FacileAnalysis)))

test_that("flm_def defines t-test and anova models on FacileBiocDataStore", {
  for (bclass in names(FBIOC)) {
    f <- FBIOC[[bclass]]
    des.ttest <- flm_def(f, "sample_type", "tumor", "normal", batch = "sex")
    expect_is(des.ttest, "FacileTtestModelDefinition",
              info = sprintf("ttest (%s)", bclass))

    des.anova <- flm_def(f, covariate = "stage", batch = "sex")
    expect_is(des.anova, "FacileAnovaModelDefinition",
              info = sprintf("anova (%s)", bclass))
  }
})

test_that("biocbox can be constructed from a FacileLinearModelDefinition", {
  for (bclass in names(FBIOC)) {
    f <- FBIOC[[bclass]]
    des.ttest <- flm_def(f, covariate = "sample_type", "tumor", "normal",
                         batch = "sex")
    des <- design(des.ttest)
    expect_equal(
      rownames(des),
      with(samples(f), paste(dataset, sample_id, sep = "__")),
      info = sprintf("rownames(design) (%s)", bclass))

    bb <- biocbox(des.ttest, filter = FALSE)
    expect_is(bb, "EList", info = bclass)
    checkmate::expect_matrix(bb[["E"]], nrows = nrow(f), ncols = ncol(f),
                             info = bclass)
    if (bclass != "EList") {
      # It would have been voomed
      checkmate::expect_matrix(bb[["weights"]], nrows = nrow(f), ncols = ncol(f),
                               info = bclass)
    }
  }
})


test_that("flm_def defines models on facile_frame from a FacileBiocDataStore", {
  for (bclass in names(FBIOC)) {
    f <- FBIOC[[bclass]]
    blca <- filter_samples(f, indication == "BLCA")
    des.ttest <- flm_def(blca, covariate = "sample_type", "tumor", "normal",
                         batch = "sex")
    expect_is(des.ttest, "FacileTtestModelDefinition", info = bclass)
    des.anova <- flm_def(blca, covariate = "stage", batch = "sex")
    expect_is(des.anova, "FacileAnovaModelDefinition", info = "bclass")

    # The `blca` sample didn't have any covariates, but thyy get spanked on
    # within flm_def so we can keep track of what the covariates were used when
    # tested. Retrieving those covariates back to blca should make it equivalent
    # again
    tsamples <- with_sample_covariates(blca, c("sample_type", "sex"))
    asamples <- with_sample_covariates(blca, c("sex", "stage"))
    expect_equal(samples(des.ttest), tsamples)
    expect_equal(samples(des.anova), asamples)
  }
})

test_that("fdge/voom works", {
  vm.facile <- FDS %>%
    filter_samples(indication == "BLCA") %>%
    flm_def("sample_type", "tumor", "normal", batch = "sex") %>%
    fdge(method = "voom")

  stats.facile <- vm.facile %>%
    tidy() %>%
    select(feature_id, name, logFC, t, pval, padj, CI.L, CI.R, B) %>%
    arrange(feature_id)

  for (bclass in .rnaseq.class) {
    f <- FBIOC[[bclass]]
    vm.bioc <- f %>%
      filter_samples(indication == "BLCA") %>%
      flm_def("sample_type", "tumor", "normal", batch = "sex") %>%
      fdge(method = "voom")
    stats.bioc <- vm.bioc %>%
      tidy() %>%
      select(colnames(stats.facile)) %>%
      arrange(feature_id)
    expect_equal(stats.bioc, stats.facile, info = bclass)
  }
})
