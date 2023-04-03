context("FacileData API: biocbox")

.classes <- c("DESeqDataSet", "DGEList", "EList", "ExpressionSet",
              "SummarizedExperiment")
.rnaseq.class <- setdiff(.classes, "EList")

if (!exists("FDS")) FDS <- FacileData::exampleFacileDataSet()
if (!exists("Y")) Y <- example_bioc_data("DGEList", efds = FDS)

BIOC <- sapply(.classes, example_bioc_data, Y = Y, simplify = FALSE)

test_that("biocbox.FacileBiocDataStore can give itself back", {
  some.samples <- collect(samples(FDS)) |> sample_frac(0.5)
  some.features <- collect(features(FDS)) |> sample_frac(0.5)

  for (bclass in names(BIOC)) {
    obj <- BIOC[[bclass]]
    f <- facilitate(obj)
    bb <- biocbox(f)
    expect_is(bb, is(obj), info = paste(bclass, "(full)"))
    expect_equal(nrow(bb), nrow(obj))
    expect_equal(ncol(bb), ncol(obj))

    bb.some <- biocbox(f, features = some.features, samples = some.samples)
    # features subset works and returned in same order as requested
    expect_equal(rownames(bb.some), some.features[["feature_id"]],
                 info = paste(bclass, "subset features"))
    # sample subset works
    expect_equal(
      with(some.samples, paste0(dataset, sample_id)),
      with(samples(bb.some), paste0(dataset, sample_id)),
      info = paste(bclass, "subset samples"))
  }
})
