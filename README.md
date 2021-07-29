
<!-- README.md is generated from README.Rmd. Please edit that file -->

# supercells

<!-- badges: start -->
[![R-CMD-check](https://github.com/Nowosad/supercells/workflows/pkgdown/badge.svg)](https://github.com/Nowosad/supercells/actions)
[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![Codecov test
coverage](https://codecov.io/gh/Nowosad/supercells/branch/master/graph/badge.svg)](https://codecov.io/gh/Nowosad/supercells?branch=master)
<!-- badges: end -->

The goal of **supercells** is to utilize the concept of superpixels to a
variety of spatial data. This package works on spatial data with one
variable (e.g., continuous raster), many variables (e.g., RGB rasters),
and spatial patterns (e.g., areas in categorical rasters). It is based
on the SLIC algorithm (Achanta et al. (2012),
<doi:10.1109/TPAMI.2012.120>), and readapts it to work with arbitrary
dissimilarity measures.

## Installation

<!-- You can install the released version of supercells from [CRAN](https://CRAN.R-project.org) with: -->

<!-- ``` r -->

<!-- install.packages("supercells") -->

<!-- ``` -->

You can install the development version from
[GitHub](https://github.com/) with:

``` r
# install.packages("remotes")
remotes::install_github("Nowosad/supercells")
```

## Example

``` r
library(supercells)
library(terra)
#> terra version 1.3.17
library(sf)
#> Linking to GEOS 3.9.0, GDAL 3.2.2, PROJ 7.2.1
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

## Documentation

See the package’s vignettes:

1.  [Superpixels of a single raster
    layer](https://nowosad.github.io/supercells/articles/articles/one_var.html)
2.  [Superpixels of an RGB
    raster](https://nowosad.github.io/supercells/articles/articles/rgb_vars.html)
    <!-- 3. [Superpixels of spatial categorical patterns](https://nowosad.github.io/supercells/articles/motifels.html) -->

## Contribution

Contributions to this package are welcome - let me know if you need
other distance measures or transformations, have any suggestions, or
spotted a bug. The preferred method of contribution is through a GitHub
pull request. Feel also free to contact us by creating [an
issue](https://github.com/nowosad/supercells/issues).
