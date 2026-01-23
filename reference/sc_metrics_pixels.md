# Pixel-level supercells metrics

Computes per-pixel spatial, value, and combined distance diagnostics

## Usage

``` r
sc_metrics_pixels(raster, x, dist_fun = "euclidean", compactness, step)
```

## Arguments

- raster:

  The input SpatRaster used to create `x`

- x:

  An sf object returned by
  [`sc_slic()`](https://jakubnowosad.com/supercells/reference/sc_slic.md)

- dist_fun:

  A distance function name or function, as in
  [`sc_slic()`](https://jakubnowosad.com/supercells/reference/sc_slic.md)

- compactness:

  A compactness value used for the supercells If missing, uses
  `attr(x, "compactness")` when available

- step:

  A step value used for the supercells If missing, uses
  `attr(x, "step")` when available

## Value

A SpatRaster with three layers:

- spatial:

  Spatial distance from each pixel to its supercell center in grid-cell
  units (row/column index distance).

- value:

  Value distance from each pixel to its supercell center in the raster
  value space.

- combined:

  Combined distance using `compactness` and `step` to scale value and
  spatial distances.

## Details

If `x` lacks `supercells`, `x`, or `y` columns, they are derived from
geometry and row order, which may differ from the original centers

## See also

[`sc_slic()`](https://jakubnowosad.com/supercells/reference/sc_slic.md),
[`sc_metrics_clusters()`](https://jakubnowosad.com/supercells/reference/sc_metrics_clusters.md),
[`sc_metrics_global()`](https://jakubnowosad.com/supercells/reference/sc_metrics_global.md)

## Examples

``` r
library(supercells)
vol = terra::rast(system.file("raster/volcano.tif", package = "supercells"))
vol_sc = sc_slic(vol, step = 8, compactness = 1)
metrics = sc_metrics_pixels(vol, vol_sc)
terra::panel(metrics, nr = 1)
```
