
<!-- README.md is generated from README.Rmd. Please edit that file -->

# supercell

<!-- badges: start -->
[![R-CMD-check](https://github.com/Nowosad/supercell/workflows/R-CMD-check/badge.svg)](https://github.com/Nowosad/supercell/actions)
[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
<!-- badges: end -->

The goal of **supercell** is to utilize the concept of superpixels to a
variety of spatial data.

## Installation

<!-- You can install the released version of supercell from [CRAN](https://CRAN.R-project.org) with: -->
<!-- ``` r -->
<!-- install.packages("supercell") -->
<!-- ``` -->

You can install the development version from
[GitHub](https://github.com/) with:

``` r
# install.packages("remotes")
remotes::install_github("Nowosad/spDataLarge")
remotes::install_github("Nowosad/supercell")
```

## Example

``` r
library(supercell)
library(terra)
#> terra version 1.1.17
srtm = rast(system.file("raster/srtm.tif", package = "spDataLarge"))
plot(srtm)
```

<img src="man/figures/README-unnamed-chunk-2-1.png" width="100%" />

``` r
srtm_slic = supercell(srtm, 23, 100, "euclidean")
plot(srtm)
plot(srtm_slic, add = TRUE, col = NA)
```

<img src="man/figures/README-unnamed-chunk-3-1.png" width="100%" />
