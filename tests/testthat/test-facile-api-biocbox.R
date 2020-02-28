context("Facile API: biocbox")

.classes <- c("DESeqDataSet", "DGEList", "EList", "ExpressionSet",
              "SummarizedExperiment")
.rnaseq.class <- setdiff(.classes, "EList")

if (!exists("FDS")) FDS <- FacileData::exampleFacileDataSet()
if (!exists("Y")) Y <- example_bioc_data("DGEList", efds = FDS)

BIOC <- sapply(.classes, example_bioc_data, Y = Y, simplify = FALSE)

test_that("biocbox returns self", {
  for (bclass in names(BIOC)) {
    obj <- BIOC[[bclass]]
    f <- facilitate(obj)
    bbox <- biocbox(f)
    expect_is(bbox, class(obj), info = bclass)
    expect_equal(colnames(bbox), colnames(obj))
    expect_equal(rownames(bbox), rownames(obj))
  }
})

test_that("biocbox with feature and sample subsets returns slimmer object", {
  sub.samples <- samples(FDS) %>% collect() %>% sample_frac(0.5)
  sub.features <- features(FDS) %>% sample_frac(0.5)
  for (bclass in names(BIOC)) {
    obj <- BIOC[[bclass]]
    f <- facilitate(obj)
    bbox <- biocbox(f, features = sub.features, samples = sub.samples)
    expect_equal(colnames(bbox), sub.samples[["sample_id"]], info = bclass)
    expect_set_equal(rownames(bbox), sub.features[["feature_id"]], info = bclass)
    for (aname in assay_names(bbox)) {
      asi <- assay_sample_info(bbox, aname)
      expect_set_equal(
        paste(asi[["dataset"]], asi[["sample_id"]]),
        paste(sub.samples[["dataset"]], sub.samples[["sample_id"]]),
        info = sprintf("assay_sample_info: %s (%s)", aname, bclass))
    }
  }
})
