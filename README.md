
<!-- README.md is generated from README.Rmd. Please edit that file -->

# FacileBioc

<!-- 
badges: start
checkout https://lazappi.github.io/clustree/ package for some badge-inspiration

[![Project Status](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active)
[![Lifecycle: Experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://www.tidyverse.org/lifecycle/#experimental)
[![Travis build status](https://travis-ci.org/facileverse/FacileBioc.svg?branch=master)](https://travis-ci.org/facileverse/FacileBioc)
[![Codecov test coverage](https://codecov.io/gh/facileverse/FacileBioc/branch/master/graph/badge.svg)](https://codecov.io/gh/facileverse/FacileBioc?branch=master)

badges: end -->

**NOTE: Development on this package has not yet begun**

The `FacileBioc` package will enable Bioconductor-standard assay
containers, like a `SummarizedExperiment`, `MultiAssayExperiment`,
`DGElist`, etc. to be used as a FacileDataStore by implementing the
FacileData API over them.

This will enable them to be used as “first-class” data-providers within
the facile ecosystem, enabling them to take advantage of all the
goodness we have on offer here.

## Example Usage

The current plan is to defined `FacileData::facilitate()` functions over
the various Bioconductor assay containers so that they return a
“facile-wrapped” version of the data container.

We would then be able to perform a tumor-vs-normal differential
expression analysis with the
[FacileAnalysis](https://github.com/facileverse/FacileAnalysis) package,
and interact with the result via `shine()`, like so:

``` r
library(FacileAnalysis)
library(edgeR)
y <- DGEList(counts = count.matrix, genes = gene.info, samples = pheno.info)
yf <- FacileBioc::facilitate(y)
dge.facile <- yf %>% 
  flm_def(group = "sample_type", numer = "tumor", denom = "normal") %>% 
  fdge(method = "voom")
shine(dge)
```

We should also still be able to use the FacileDGEList object `yf` as a
normal `DGEList` as well, ie.

``` r
yf <- calcNormFactors(yf)
yf <- estimateDisp(yf, model.matrix(~ sample_type, yf$samples), robust = TRUE)
fit <- glmQLFit(yf, yf$design, robust = TRUE)
dge.standard <- glmQLFTest(fit, coef = "sample_typetumor") 
```

… or something like that …
