# Plot iteration diagnostics

Plot mean distance across iterations for a supercells run

## Usage

``` r
sc_plot_iter_diagnostics(x)
```

## Arguments

- x:

  A supercells object with an `iter_diagnostics` attribute, or a
  diagnostics list containing `mean_distance`

## Value

Invisibly returns `TRUE` when a plot is created

## See also

[`sc_slic()`](https://jakubnowosad.com/supercells/reference/sc_slic.md),
[`sc_slic_points()`](https://jakubnowosad.com/supercells/reference/sc_slic_points.md),
[`sc_slic_raster()`](https://jakubnowosad.com/supercells/reference/sc_slic_raster.md)

## Examples

``` r
library(supercells)
vol = terra::rast(system.file("raster/volcano.tif", package = "supercells"))
vol_sc = sc_slic(vol, step = 8, compactness = 1, iter_diagnostics = TRUE)
sc_plot_iter_diagnostics(vol_sc)
```
