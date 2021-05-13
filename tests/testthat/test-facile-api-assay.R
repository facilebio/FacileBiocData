context("FacileData API: assay-level query and retrieval")


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
    f <- facilitate(obj, run_vst = FALSE)
    expected_assays <- assay_names(obj)
    if (is(obj, "DESeqDataSet")) {
      expected_assays <- c(expected_assays, "normcounts")
    }
    checkmate::expect_set_equal(assay_names(f), expected_assays, info = bclass)
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
    f <- facilitate(obj, run_vst = FALSE)
    ainfo <- assay_info(f)
    expect_s3_class(ainfo, "data.frame")

    # check each column is of expected data type
    for (cname in names(expected.cols)) {
      ctype <- expected.cols[cname]
      info <- sprintf("column '%s' is not type '%s' from container '%s'",
                      cname, ctype, bclass)
      expect_is(ainfo[[cname]], ctype, info = info)
    }

    expect_equal(ainfo[["feature_type"]][1L], "entrez", info = bclass)
    expect_equal(ainfo[["nfeatures"]][1L], unname(nrow(obj)), info = bclass)
    checkmate::expect_set_equal(
      ainfo[["assay"]],
      assay_names(f),
      info = bclass)
  }
})

test_that("(fetch|with)_assay_data retrieval works across containers", {
  # we need to sample the same genes all the time because the tolerance
  # setting required for DESeq2 normalized counts to match edgeR:cpm counts
  # may jump around for different sets of genes.
  set.seed(42)
  # This test is restricted to rnaseq containers for now
  features.all <- features(FDS)
  features.some <- dplyr::sample_n(features.all, 5)
  samples.all <- samples(FDS) %>% collect()
  samples.some <- dplyr::sample_n(samples.all, 10)

  # The names of the assay will differ accross bioc data container types,
  # so we remove that column from these results
  fds.res <- list(
    tidy.all = FDS %>%
      fetch_assay_data(features.some, samples.all) %>%
      select(-assay) %>%
      arrange(sample_id, feature_id),
    tidy.some = FDS %>%
      fetch_assay_data(features.some, samples.some) %>%
      select(-assay) %>%
      arrange(sample_id, feature_id),
    tidy.some.fids = FDS %>%
      fetch_assay_data(features.some$feature_id, samples.some) %>%
      select(-assay) %>%
      arrange(sample_id, feature_id),
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

  # See the next test for SummarizedExperiment
  for (bclass in setdiff(.rnaseq.class, "SummarizedExperiment")) {
    obj <- BIOC[[bclass]]

    f <- facilitate(obj, assay_type = "rnaseq", run_vst = FALSE)
    bsamples.all <- samples(f)
    bsamples.some <- semi_join(bsamples.all, samples.some,
                               by = c("dataset", "sample_id"))
    normalized <- TRUE

    bioc.res <- list(
      # exclude the assay name from the tidy'd results because they will differ
      # across containers
      tidy.all = f %>%
        fetch_assay_data(features.some, bsamples.all) %>%
        select(-assay) %>%
        arrange(sample_id, feature_id),
      tidy.some = f %>%
        fetch_assay_data(features.some, bsamples.some) %>%
        select(-assay) %>%
        arrange(sample_id, feature_id),
      tidy.some.fids = f %>%
        fetch_assay_data(features.some$feature_id, bsamples.some) %>%
        select(-assay) %>%
        arrange(sample_id, feature_id),
      tidy.with = bsamples.some %>%
        with_assay_data(features.some, normalized = normalized) %>%
        arrange(sample_id),
      matrix.all = f %>%
        fetch_assay_data(features.some, bsamples.all, as.matrix = TRUE),
      matrix.some.norm = f %>%
        fetch_assay_data(features.some, bsamples.some, normalized = normalized,
                         as.matrix = TRUE))

    for (comp in names(bioc.res)) {
      bres <- bioc.res[[comp]]
      fres <- fds.res[[comp]]
      is.tidy <- grepl("tidy\\.", comp)
      info <- sprintf("[%s] %s (%s)", bclass, sub("^.*?\\.", "", comp),
                      sub("\\..*$", "", comp))

      expect_is(fds(bres), is(f), info = info) # Ensure results were generated
      expect_is(fds(bres), is(obj))            # from correct container type.
      expect_is(bres, is(fres), info = info)   # Results are the same type.

      # normalization from DESeqDataSet uses DESeq2::count(x, normalized = TRUE)
      # which is different than the cpms that we are testing against
      if (bclass == "DESeqDataSet") {
        tolerance <- 2.22
      } else {
        tolerance <- testthat::testthat_tolerance()
      }
      if (is.tidy) {
        # expect_equal(bres, fres, info = info)
        # if this is from deseqdataset, normfactor and libsize came along
        # for the ride
        ecols <- colnames(fres)
        checkmate::expect_subset(ecols, colnames(bres), info = info)
        # For some reason the default expect_equal.tbl_df would randomly break
        # on SummarizedExperiment results. Manually expect_equal()-ing the
        # columns of the tbl's against each other always returned TRUE, though,
        # and I think it boils down to this:
        # https://github.com/tidyverse/dplyr/issues/2751
        # expect_equal.tbl_df has been removed from dplyr 1.0, but we needed
        # this to work before that, so .....................
        expect_equal(
          as.data.frame(bres[, ecols]),
          as.data.frame(fres),
          tolerance = tolerance,
          info = info,
          check.attributes = FALSE)
      } else {
        checkmate::expect_set_equal(colnames(bres), colnames(fres))
        checkmate::expect_set_equal(rownames(bres), rownames(fres))
        bres <- bres[rownames(fres), colnames(fres)]
        expect_equal(
          as.data.frame(bres),
          as.data.frame(fres),
          tolerance = tolerance,
          check.attributes = FALSE,
          info = info)
      }
    }
  }
})

# test_that("Stress test facile api against SummarizedExperiment", {
#   # Tests against SummarizedExperiment seems to randomly fail on line 150
#   # in the "tidy" tests ... need to loop on this to find it. The other
#   # containers never seem to fail,
#
#   obj <- BIOC[["SummarizedExperiment"]]
#   f <- facilitate(obj, assay_type = "rnaseq")
#   normalized <- TRUE
#
#   features.all <- features(FDS)
#   samples.all <- samples(FDS) %>% collect()
#
#   for (i in 1:10) {
#     features.some <- dplyr::sample_n(features.all, 5)
#     samples.some <- dplyr::sample_n(samples.all, 10)
#
#     # The names of the assay will differ accross bioc data container types,
#     # so we remove that column from these results
#     fds.res <- list(
#       tidy.all = FDS %>%
#         fetch_assay_data(features.some, samples.all) %>%
#         select(-assay) %>%
#         arrange(sample_id, feature_id),
#       tidy.some = FDS %>%
#         fetch_assay_data(features.some, samples.some) %>%
#         select(-assay) %>%
#         arrange(sample_id, feature_id),
#       tidy.some.fids = FDS %>%
#         fetch_assay_data(features.some$feature_id, samples.some) %>%
#         select(-assay) %>%
#         arrange(sample_id, feature_id),
#       # Exercising the `with_` call here simultaneously tests the with_
#       # decoration functionality as well as the normalization procedure, since
#       # the default for `with_assay_data` is `normalized = TRUE`
#       tidy.with = samples.some %>%
#         with_assay_data(features.some) %>%
#         arrange(sample_id),
#       matrix.all = FDS %>%
#         fetch_assay_data(features.some, samples.all, as.matrix = TRUE),
#       matrix.some.norm = FDS %>%
#         fetch_assay_data(features.some, samples.some, as.matrix = TRUE,
#                          normalized = TRUE))
#
#
#     bsamples.all <- samples(f)
#     bsamples.some <- semi_join(bsamples.all, samples.some,
#                                by = c("dataset", "sample_id"))
#
#     bioc.res <- list(
#       # exclude the assay name the tidied reuslts because they will differ
#       # across containers
#       tidy.all = f %>%
#         fetch_assay_data(features.some, bsamples.all) %>%
#         select(-assay) %>%
#         arrange(sample_id, feature_id),
#       tidy.some = f %>%
#         fetch_assay_data(features.some, bsamples.some) %>%
#         select(-assay) %>%
#         arrange(sample_id, feature_id),
#       tidy.some.fids = f %>%
#         fetch_assay_data(features.some$feature_id, bsamples.some) %>%
#         select(-assay) %>%
#         arrange(sample_id, feature_id),
#       tidy.with = bsamples.some %>%
#         with_assay_data(features.some, normalized = normalized) %>%
#         arrange(sample_id),
#       matrix.all = f %>%
#         fetch_assay_data(features.some, bsamples.all, as.matrix = TRUE),
#       matrix.some.norm = f %>%
#         fetch_assay_data(features.some, bsamples.some, normalized = normalized,
#                          as.matrix = TRUE))
#
#     for (comp in names(bioc.res)) {
#       bres <- bioc.res[[comp]]
#       fres <- fds.res[[comp]]
#       is.tidy <- grepl("tidy\\.", comp)
#       info <- sprintf("[%s] %s (%s)", bclass, sub("^.*?\\.", "", comp),
#                       sub("\\..*$", "", comp))
#
#       expect_is(fds(bres), is(f), info = info) # Ensure results were generated
#       expect_is(fds(bres), is(obj))            # from correct container type.
#       expect_is(bres, is(fres), info = info)   # Results are the same type.
#       if (is.tidy) {
#         ecols <- colnames(fres)
#         checkmate::expect_subset(ecols, colnames(bres), info = info)
#         expect_equal(
#           as.data.frame(bres[, ecols]),
#           as.data.frame(fres),
#           info = info,
#           check.attributes = FALSE)
#       } else {
#         checkmate::expect_set_equal(colnames(bres), colnames(fres))
#         checkmate::expect_set_equal(rownames(bres), rownames(fres))
#         bres <- bres[rownames(fres), colnames(fres)]
#         expect_equal(bres, fres, check.attributes = FALSE, info = info)
#       }
#     }
#   }
#
# })
