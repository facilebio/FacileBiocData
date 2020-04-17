library("testthat")
library("checkmate")
library("dplyr")
library("magrittr")
library("FacileBiocData")

FDS <- FacileData::exampleFacileDataSet()
Y <- example_bioc_data("DGEList", FDS)

test_check("FacileBiocData")
