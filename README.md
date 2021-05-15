
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
#> terra version 1.2.10
library(sf)
#> Linking to GEOS 3.8.1, GDAL 3.1.4, PROJ 6.3.2
vol = rast(system.file("raster/volcano.tif", package = "supercells"))
plot(vol)
```

<img src="man/figures/README-unnamed-chunk-2-1.png" width="100%" />

``` r
vol_slic1 = supercells(vol, k = 50, compactness = 1)
plot(vol)
plot(st_geometry(vol_slic1), add = TRUE, lwd = 0.2)
```

<img src="man/figures/README-unnamed-chunk-3-1.png" width="100%" />

``` r
vol_slic2 = supercells(vol, k = 50, compactness = 1, dist_fun = "euclidean", 
                        clean = FALSE, iter = 10)
plot(vol)
plot(st_geometry(vol_slic2), add = TRUE, lwd = 0.2)
```

<img src="man/figures/README-unnamed-chunk-4-1.png" width="100%" />

## Contribution

Contributions to this package are welcome - let me know if you need
other distance measures or transformations, have any suggestions, or
spotted a bug. The preferred method of contribution is through a GitHub
pull request. Feel also free to contact us by creating [an
issue](https://github.com/nowosad/supercells/issues).