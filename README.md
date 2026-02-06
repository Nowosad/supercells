
<!-- README.md is generated from README.Rmd. Please edit that file -->

# supercells <img src="man/figures/logo.png" align="right" width="150" />

<!-- badges: start -->

[![CRAN
status](https://www.r-pkg.org/badges/version/supercells)](https://CRAN.R-project.org/package=supercells)
[![R-CMD-check](https://github.com/Nowosad/supercells/workflows/pkgdown/badge.svg)](https://github.com/Nowosad/supercells/actions)
[![Codecov test
coverage](https://codecov.io/gh/Nowosad/supercells/branch/master/graph/badge.svg)](https://app.codecov.io/gh/Nowosad/supercells?branch=master)
[![CRAN RStudio mirror
downloads](http://cranlogs.r-pkg.org/badges/supercells)](https://cran.r-project.org/package=supercells)
<!-- badges: end -->

The goal of **supercells** is to bring the idea of superpixels into
spatial analysis - grouping raster cells (not just image pixels) into
meaningful, coherent regions. This package works on spatial data with
one variable (e.g., continuous raster), many variables (e.g., RGB
rasters), and spatial patterns (e.g., areas in categorical rasters). It
is based on the SLIC algorithm (Achanta et al. (2012)), adapted to work
with arbitrary dissimilarity measures.

## Installation

You can install the released version of supercells from
[CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("supercells")
```

You can install the development version from
[GitHub](https://github.com/) with:

``` r
install.packages("supercells", repos = "https://nowosad.r-universe.dev")
```

<!-- ``` r -->
<!-- # install.packages("remotes") -->
<!-- remotes::install_github("Nowosad/supercells") -->
<!-- ``` -->

## Example

``` r
library(supercells)
library(terra)
#> terra 1.6.49
library(sf)
#> Linking to GEOS 3.11.0, GDAL 3.5.2, PROJ 9.0.1; sf_use_s2() is TRUE
vol = rast(system.file("raster/volcano.tif", package = "supercells"))
plot(vol)
```

<img src="man/figures/README-unnamed-chunk-2-1.png" width="100%" />

``` r
vol_slic1 = sc_slic(vol, step = 8, compactness = 1)
plot(vol)
plot(st_geometry(vol_slic1), add = TRUE, lwd = 0.5, border = "red")
```

<img src="man/figures/README-unnamed-chunk-3-1.png" width="100%" />

## Documentation

The upcoming version 2 (development) introduces a refreshed API with
`sc_`-prefixed functions (for example, `sc_slic()`), plus new tooling for
parameter tuning and evaluation. Start with these v2 vignettes:

1.  [Main changes since version
    1.0.0](https://jakubnowosad.com/supercells/articles/v2-changes-since-v1.html)
2.  [Introduction to
    supercells](https://jakubnowosad.com/supercells/articles/v2-intro.html)
3.  [Choosing parameters for
    supercells](https://jakubnowosad.com/supercells/articles/v2-parameters.html)
4.  [Evaluation and
    diagnostics](https://jakubnowosad.com/supercells/articles/v2-evaluation.html)

Watch the presentations about this package and some related ideas:

1.  *Spatial segmentation in R using the supercells package*,
    2021-09-02, OpenGeoHub Summer School -
    [slides](https://jakubnowosad.com/ogh2021/),
    [video](https://doi.org/10.5446/54880)
2.  *Generalizing the Simple Linear Iterative Clustering (SLIC)
    superpixels*, 2021-09-28, GIScience 2021 -
    [slides](https://jakubnowosad.com/giscience-2021/),
    [video](https://www.youtube.com/watch?v=AlyEFkyKLUw&t=2018s)

Read the related article:

1.  Nowosad, J., Stepinski, T. (2022). Extended SLIC superpixels
    algorithm for applications to non-imagery geospatial rasters.
    International Journal of Applied Earth Observation and
    Geoinformation, <https://doi.org/10.1016/j.jag.2022.102935>

## Contribution

Contributions to this package are welcome - let us know if you need
other distance measures or transformations, have any suggestions, or
spotted a bug. The preferred method of contribution is through a GitHub
pull request. Feel also free to contact us by creating [an
issue](https://github.com/nowosad/supercells/issues).
