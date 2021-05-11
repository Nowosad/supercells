
<!-- README.md is generated from README.Rmd. Please edit that file -->

# supercells

<!-- badges: start -->
[![R-CMD-check](https://github.com/Nowosad/supercells/workflows/R-CMD-check/badge.svg)](https://github.com/Nowosad/supercells/actions)
[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![Codecov test
coverage](https://codecov.io/gh/Nowosad/supercells/branch/master/graph/badge.svg)](https://codecov.io/gh/Nowosad/supercells?branch=master)
<!-- badges: end -->

The goal of **supercells** is to utilize the concept of superpixels to a
variety of spatial data.

## Installation

<!-- You can install the released version of supercells from [CRAN](https://CRAN.R-project.org) with: -->
<!-- ``` r -->
<!-- install.packages("supercells") -->
<!-- ``` -->

You can install the development version from
[GitHub](https://github.com/) with:

``` r
# install.packages("remotes")
remotes::install_github("Nowosad/spDataLarge")
remotes::install_github("Nowosad/supercells")
```

## Example

``` r
library(supercells)
library(terra)
#> terra version 1.1.17
library(sf)
#> Linking to GEOS 3.8.1, GDAL 3.1.4, PROJ 6.3.2
srtm = rast(system.file("raster/srtm.tif", package = "spDataLarge"))
plot(srtm)
```

<img src="man/figures/README-unnamed-chunk-2-1.png" width="100%" />

``` r
srtm_slic1 = supercells(srtm, k = 500, compactness = 1)
plot(srtm)
plot(st_geometry(srtm_slic1), add = TRUE, lwd = 0.2)
```

<img src="man/figures/README-unnamed-chunk-3-1.png" width="100%" />

``` r
srtm_slic2 = supercells(srtm, k = 500, compactness = 1, dist_fun = "euclidean", 
                        clean = FALSE, iter = 10)
plot(srtm)
plot(st_geometry(srtm_slic2), add = TRUE, lwd = 0.2)
```

<img src="man/figures/README-unnamed-chunk-4-1.png" width="100%" />
