context("FacileAnalysis: PCA")

.classes <- c("DESeqDataSet", "DGEList", "EList", "ExpressionSet",
              "SummarizedExperiment")

if (!exists("FDS")) FDS <- FacileData::exampleFacileDataSet()
if (!exists("Y")) Y <- example_bioc_data("DGEList", efds = FDS)

BIOC <- sapply(.classes, example_bioc_data, Y = Y, simplify = FALSE)
FBIOC <- lapply(BIOC, function(b) {
  assay_type <- if (is(b, "EList")) "lognorm" else "rnaseq"
  facilitate(b, assay_type = assay_type, run_vst = FALSE)
})

suppressPackageStartupMessages(suppressWarnings(library(FacileAnalysis)))

test_that("fpca works", {
  set.seed(123)
  # prior.count set to match what is used in voom, since we are testing
  # an EList, which was generated within this package by voom
  pca.facile <- FDS |>
    filter_samples(indication == "BLCA") |>
    FacileAnalysis::fpca(ntop = 1000, prior.count = 0.5)

  sig.facile <- FacileAnalysis::signature(pca.facile) |> tidy()
  top.facile <- sig.facile |>
    filter(weight >= quantile(sig.facile$weight, 0.75)) |>
    distinct(feature_id)
  for (bclass in names(FBIOC)) {
    f <- FBIOC[[bclass]]
    set.seed(123)
    # There is an issue casting the survival columns we've got in the
    # dataset, let's ignore that for now
    pca.bioc <- expect_warning({
      f |>
        filter_samples(indication == "BLCA") |>
        FacileAnalysis::fpca(features = features(pca.facile), prior.count = 0.5)
    }, "PFS.*conversion")

    # percent variance explained should be roughly similar
    expect_equal(pca.bioc$percent_var, pca.facile$percent_var,
                 tolerance = 0.01,
                 info = sprintf("var explained (%s)", bclass))

    # check top loading genes are about the same, ie. the top loaded genes
    # from the bioc container covers > 95% of the top genes from the FDS
    sig.bioc <- FacileAnalysis::signature(pca.bioc) |> tidy()
    f.recovered <- intersect(sig.bioc$feature_id, top.facile$feature_id)
    f.fraction <- mean(top.facile$feature_id %in% f.recovered)
    expect_true(f.fraction >= 0.98,
                info = sprintf("fraction top genes recovered (%s)", bclass))
  }
})
