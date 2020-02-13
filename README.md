
<!-- README.md is generated from README.Rmd. Please edit that file -->

# FacileBiocData

<!-- 
badges: start
checkout https://lazappi.github.io/clustree/ package for some badge-inspiration

[![Project Status](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active)
[![Lifecycle: Experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://www.tidyverse.org/lifecycle/#experimental)
[![Travis build status](https://travis-ci.org/facilebio/FacileBiocData.svg?branch=master)](https://travis-ci.org/facilebio/FacileBiocData)
[![Codecov test coverage](https://codecov.io/gh/facilebio/FacileBioc/branch/master/graph/badge.svg)](https://codecov.io/gh/facilebio/FacileBiocData?branch=master)

badges: end -->

**NOTE**: This package is incomplete. An alpha version should be
available for use in early March, 2020.

The `FacileBiocData` package will define the FacileData API over
Bioconductor-standard assay containers, like a `SummarizedExperiment`,
`MultiAssayExperiment`, `DGEList`, etc.

This will enable them to be used as “first-class” data-providers within
the facile ecosystem, so that they can take advantage of all the
goodness we have on offer here.

## Example Usage

The current plan is to implement `FacileData::facilitate()` functions
over the the assay containers (eg. `facilitate.DGEList()`) so that they
return a “facile-wrapped” version of the data container.

``` r
library(FacileBiocData)

efds <- FacileData::exampleFacileDataSet()
y <- FacileData::as.DGEList(efds)
yf <- facilitate(y)
```

Once we’ve decorated the base Bioconductor assay container, we can use
it within the facile.bio ecosystem. For instance, we can now perform a
differential gene expression analysis using the \`FacileAnalysis
package.

``` r
library(FacileAnalysis)
dge.facile <- yf %>% 
  flm_def(group = "sample_type", numer = "tumor", denom = "normal") %>% 
  fdge(method = "edgeR-qlf")
```

``` r
shine(dge.facile)
```

We should still be able to use the `FacileDGEList` object (`yf`) as a
normal `DGEList`, so that the “normal” edgeR moves would work as well.

``` r
genes <- features(dge.facile)$feature_id # restrict to genes used in fdge()
yf <- edgeR::calcNormFactors(yf[genes,,keep.lib.sizes = FALSE])
yf <- edgeR::estimateDisp(yf, model.matrix(~ sample_type, data = yf$samples))
fit <- edgeR::glmQLFit(yf, yf$design, robust = TRUE)
results <- edgeR::glmQLFTest(fit, coef = "sample_typetumor")
dge.standard <- edgeR::topTags(results, n = Inf, sort.by = "none")
```

And these two results are equivalent

``` r
all.equal(tidy(dge.facile)$logFC, dge.standard$logFC)
```
