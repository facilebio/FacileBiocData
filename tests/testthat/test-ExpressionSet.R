context("ExpressionSet")

library(Biobase)

test_that("facilitate.ExpressionSet works", {
  efds <- FacileData::exampleFacileDataSet()
  y <- FacileData::as.DGEList(efds)
  eset <- ExpressionSet(y$counts)
  pData(eset) <- y$samples
  fData(eset) <- y$genes

  feset <- facilitate(eset)
})
