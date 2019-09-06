
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

The `FacileBioc` package will define the FacileData API over
Bioconductor-standard assay containers, like a `SummarizedExperiment`,
`MultiAssayExperiment`, `DGEList`, etc.

This will enable them to be used as “first-class” data-providers within
the facile ecosystem, so that they can take advantage of all the
goodness we have on offer here.

## Example Usage

The current plan is to implement `FacileData::facilitate()` functions
over the the assay containers (eg. `facilitate.DGEList()`) so that they
return a “facile-wrapped” version of the data container.

Imagine we had a `DGEList` object (`y`), with RNA-seq data from an
experiment with “normal” and “tumor” samples. A differential expression
analysis contrasting the two sample types using the
[FacileAnalysis](https://github.com/facileverse/FacileAnalysis) package
would like like so:

``` r
library(FacileAnalysis)
library(edgeR)
y <- DGEList(counts = count.matrix, genes = gene.info, samples = sample.info)
yf <- FacileBioc::facilitate(y)
dge.facile <- yf %>% 
  flm_def(group = "sample_type", numer = "tumor", denom = "normal") %>% 
  fdge(method = "voom")
shine(dge.facile)
```

We should still be able to use the `FacileDGEList` object (`yf`) as a
normal `DGEList`, so that the “normal” edgeR moves would work as well:

``` r
yf <- calcNormFactors(yf)
yf <- estimateDisp(yf, model.matrix(~ sample_type, yf$samples), robust = TRUE)
fit <- glmQLFit(yf, yf$design, robust = TRUE)
dge.standard <- glmQLFTest(fit, coef = "sample_typetumor") 
```

… or something like that …
