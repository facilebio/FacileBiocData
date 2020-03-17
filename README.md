
<!-- README.md is generated from README.Rmd. Please edit that file -->

# FacileBiocData

<!-- 
badges: start
checkout https://lazappi.github.io/clustree/ package for some badge-inspiration

[![Travis build status](https://travis-ci.org/facilebio/FacileBiocData.svg?branch=master)](https://travis-ci.org/facilebio/FacileBiocData)
[![Codecov test coverage](https://codecov.io/gh/facilebio/FacileBioc/branch/master/graph/badge.svg)](https://codecov.io/gh/facilebio/FacileBiocData?branch=master)

badges: end -->

[![Project
Status](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active)
[![Lifecycle:
Maturing](https://img.shields.io/badge/lifecycle-maturing-blue.svg)](https://www.tidyverse.org/lifecycle/#maturing)

The `FacileBiocData` package enables the use of Bioconductor-standard
data containers, like a `SummarizedExperiment`, `DGEList`,
`DESeqDataSet`, etc. as “first-class” data-providers within the facile
ecosystem.

## Example Usage

The user simply needs to call the `facilitate` function on their data
container in order to make its data available via the facile API, so
that it can be analyzed within the facile framework.

``` r
library(FacileBiocData)
data("airway", package = "airway")
airway.facile <- facilitate(airway, assay_type = "rnaseq")
```

We can now use `airway.facile` as a first-class data-providedr within
the facile framework. For instance, we can use the
[FacileAnalysis](https://facilebio.github.io/FacileAnalysis/) to perform
a differential expression analysis using the edgeR or limma based
framework:

``` r
library(FacileAnalysis)
dge.facile <- airway.facile %>% 
  flm_def("dex", numer = "trt", denom = "untrt", batch = "cell") %>% 
  fdge(method = "voom")
```

We can extract the statistics from the `fdge` result:

``` r
tidy(dge.facile) %>% 
  select(feature_id, logFC, pval, padj) %>% 
  arrange(pval) %>% 
  head()
#> # A tibble: 6 x 4
#>   feature_id      logFC     pval        padj
#>   <chr>           <dbl>    <dbl>       <dbl>
#> 1 ENSG00000165995  3.28 4.29e-11 0.000000684
#> 2 ENSG00000179593  8.06 3.14e-10 0.00000152 
#> 3 ENSG00000120129  2.94 5.19e-10 0.00000152 
#> 4 ENSG00000152583  4.56 5.82e-10 0.00000152 
#> 5 ENSG00000162493  1.88 6.18e-10 0.00000152 
#> 6 ENSG00000157214  1.97 7.07e-10 0.00000152
```

Produce an interactive visual (via using plotly/htmlwidgets) from one of
the results using `viz()`

``` r
viz(dge.facile, "ENSG00000165995")
```

<img src="man/figures/README-viz-fdge.png" width="50%" />

Or, finally, launch a shiny gadget over the `fdge()` result so that we
can interactively explore the differential expression result in all of
its glory:

``` r
shine(dge.facile)
```

<img src="man/figures/README-shine-fdge.png" width="75%" />

You can refer to the [RNA-seq analysis
vignette](https://facilebio.github.io/FacileAnalysis/articles/FacileAnalysis-RNAseq.html)
vignette in the
[FacileAnalysis](https://facilebio.github.io/FacileAnalysis/) package in
order to learn how you can interactively analyze and explore RNA-seq
data in the facile.bio framework.
