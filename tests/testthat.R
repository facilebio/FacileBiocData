library("testthat")
library("checkmate")
library("dplyr")
library("magrittr")
library("FacileBiocData")

FDS <- FacileData::exampleFacileDataSet()
Y <- local({
  y <- FacileData::as.DGEList(FDS)
  y[["samples"]][["samid"]] <- NULL
  colnames(y) <- y$samples$sample_id
  y[["samples"]][["group"]] <- y$samples$sample_type
  y
})

test_package("FacileAnalysis")

