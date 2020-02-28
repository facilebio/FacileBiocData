context("Facile API: assay-level query and retrieval")


.classes <- c("DESeqDataSet", "DGEList", "EList", "ExpressionSet",
              "SummarizedExperiment")
.rnaseq.class <- setdiff(.classes, "EList")

if (!exists("FDS")) FDS <- FacileData::exampleFacileDataSet()
if (!exists("Y")) Y <- example_bioc_data("DGEList", efds = FDS)

BIOC <- sapply(.classes, example_bioc_data, Y = Y, simplify = FALSE)

test_that("assay_names is defined for base Bioconductor classes", {
  for (bclass in names(BIOC)) {
    obj <- BIOC[[bclass]]
    fn.name <- paste0("assay_names.", class(obj)[1L])
    checkmate::expect_function(getFunction(fn.name), info = bclass)
    anames <- assay_names(obj)
    checkmate::expect_character(anames, min.len = 1L, info = bclass)
  }
})

test_that("assay_names returns names of assays in the facilitated container", {
  for (bclass in names(BIOC)) {
    obj <- BIOC[[bclass]]
    f <- facilitate(obj)
    checkmate::expect_set_equal(assay_names(f), assay_names(obj), info = bclass)
  }
})

test_that("assay_info returns legit metadata for all containers", {
  expected.cols <- c(
    assay = "character",
    assay_type = "character",
    feature_type = "character",
    description = "character",
    nfeatures = "integer",
    storage_mode = "character")

  for (bclass in names(BIOC)) {
    obj <- BIOC[[bclass]]
    f <- facilitate(obj)
    ainfo <- assay_info(f)
    expect_s3_class(ainfo, "data.frame")

    # check each column is of expected data type
    for (cname in names(expected.cols)) {
      ctype <- expected.cols[cname]
      info <- sprintf("column '%s' is not type '%s' from container '%s'",
                      cname, ctype, bclass)
      expect_is(ainfo[[cname]], ctype, info = info)
    }

    expect_equal(ainfo[["feature_type"]][1L], "entrez", info = blcass)
    expect_equal(ainfo[["nfeatures"]][1L], unname(nrow(obj)), info = bclass)
    checkmate::expect_set_equal(ainfo[["assay"]], assay_names(f), info = bclass)
  }
})

test_that("(fetch|with)_assay_data retrieval works across containers", {
  # This test is restricted to rnaseq containers for now
  features.all <- features(FDS)
  features.some <- sample_n(features.all, 5)
  samples.all <- samples(FDS) %>% collect()
  samples.some <- sample_n(samples.all, 10)

  # The names of the assay will differ accross bioc data container types,
  # so we remove that column from these results
  fds.res <- list(
    tidy.all = FDS %>%
      fetch_assay_data(features.some, samples.all) %>%
      select(-assay),
    tidy.some = FDS %>%
      fetch_assay_data(features.some, samples.some) %>%
      select(-assay),
    # Exercising the `with_` call here simultaneously tests the with_
    # decoration functionality as well as the normalization procedure, since
    # the default for `with_assay_data` is `normalized = TRUE`
    tidy.with = samples.some %>%
      with_assay_data(features.some) %>%
      arrange(sample_id),
    matrix.all = FDS %>%
      fetch_assay_data(features.some, samples.all, as.matrix = TRUE),
    matrix.some.norm = FDS %>%
      fetch_assay_data(features.some, samples.some, as.matrix = TRUE,
                       normalized = TRUE))

  for (bclass in .rnaseq.class) {
    obj <- BIOC[[bclass]]
    f <- facilitate(obj)
    bsamples.all <- samples(f)
    bsamples.some <- semi_join(bsamples.all, samples.some,
                               by = c("dataset", "sample_id"))

    bioc.res <- list(
      # exclude the assay name the tidied reuslts because they will differ
      # across containers
      tidy.all = f %>%
        fetch_assay_data(features.some, bsamples.all) %>%
        select(-assay),
      tidy.some = f %>%
        fetch_assay_data(features.some, bsamples.some) %>%
        select(-assay),
      tidy.with = bsamples.some %>%
        with_assay_data(features.some) %>%
        arrange(sample_id),
      matrix.all = f %>%
        fetch_assay_data(features.some, bsamples.all, as.matrix = TRUE),
      matrix.some.norm = f %>%
        fetch_assay_data(features.some, bsamples.some, as.matrix = TRUE,
                         normalized = TRUE))

    for (comp in names(bioc.res)) {
      bres <- bioc.res[[comp]]
      fres <- fds.res[[comp]]
      is.tidy <- grepl("tidy\\.", comp)
      info <- sprintf("[%s] %s (%s)", bclass, sub("^.*?\\.", "", comp),
                      sub("\\..*$", "", comp))

      expect_is(fds(bres), is(f), info = info) # Ensure results were generated
      expect_is(fds(bres), is(obj))            # from correct container type.
      expect_is(bres, is(fres), info = info)   # Results are the same type.
      if (is.tidy) {
        expect_equal(bres, fres, info = info)
      } else {
        expect_set_equal(colnames(bres), colnames(fres))
        expect_set_equal(rownames(bres), rownames(fres))
        bres <- bres[rownames(fres), colnames(fres)]
        expect_equal(bres, fres, check.attributes = FALSE, info = info)
      }
    }
  }
})
