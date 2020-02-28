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
  adat.all.fds <- FDS %>%
    fetch_assay_data(features.some, samples.all) %>%
    select(-assay)
  adat.some.fds <- FDS %>%
    fetch_assay_data(features.some, samples.some) %>%
    select(-assay)

  # Exercising the `with_` call here simultaneously tests the with_
  # decoration functionality as well as the normalization procedure, since
  # the default for `with_assay_data` is `normalized = TRUE`
  with.adat.fds <- samples.some %>%
    with_assay_data(features.some, assay_name = "rnaseq") %>%
    arrange(sample_id)

  for (bclass in .rnaseq.class) {
    obj <- BIOC[[bclass]]
    f <- facilitate(obj)
    bsamples.all <- samples(f)
    bsamples.some <- semi_join(bsamples.all, samples.some,
                               by = c("dataset", "sample_id"))

    adat.all.bioc <- f %>%
      fetch_assay_data(features.some, bsamples.all) %>%
      select(-assay)
    adat.some.bioc <- f %>%
      fetch_assay_data(features.some, bsamples.some) %>%
      select(-assay)

    # The names of the assays are different among containers, so we explicitly
    # do not test those
    expect_equal(adat.all.bioc, adat.all.fds, info = bclass)
    expect_equal(adat.some.bioc, adat.some.fds, info = bclass)

    # The order of the samples returned isn't guaranteed, so we force them
    # to be lexicographical order just so that we check the values are the
    # same
    with.adat.bioc <- bsamples.some %>%
      with_assay_data(features.some, assay_name = default_assay(f)) %>%
      arrange(sample_id)

    expect_equal(with.adat.bioc, with.adat.fds, info = bclass)

    # Let's just double-check we have been checking results from the right
    # bioconductor container
    expect_equal(class(fds(with.adat.bioc)), class(f), info = bclass)
  }
})
