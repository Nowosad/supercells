# Global supercells metrics

Computes global distance diagnostics for supercells

## Usage

``` r
sc_metrics_global(raster, x, dist_fun = "euclidean", compactness, step)
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

A data.frame with a single row of global metrics and columns:

- step:

  Step size used to generate supercells.

- compactness:

  Compactness value used to generate supercells.

- n_supercells:

  Number of supercells with at least one non-missing pixel.

- mean_value_dist:

  Mean per-supercell value distance from pixels to their supercell
  centers, averaged across supercells.

- mean_spatial_dist:

  Mean per-supercell spatial distance from pixels to their supercell
  centers, averaged across supercells; units are grid cells (row/column
  index distance), not map units.

- mean_combined_dist:

  Mean per-supercell combined distance, computed from value and spatial
  distances using `compactness` and `step`, averaged across supercells.

- compactness_ratio_mean:

  Mean ratio of scaled value distance to scaled spatial distance,
  averaged across supercells; `NA` when `compactness` or `step` is zero.

## Details

Requires `x` with metadata columns (`supercells`, `x`, `y`) If they are
missing, they are derived from geometry and row order Set
`metadata = TRUE` when calling
[`sc_slic()`](https://jakubnowosad.com/supercells/reference/sc_slic.md)
or
[`supercells()`](https://jakubnowosad.com/supercells/reference/supercells.md)

## See also

[`sc_slic()`](https://jakubnowosad.com/supercells/reference/sc_slic.md),
[`sc_metrics_pixels()`](https://jakubnowosad.com/supercells/reference/sc_metrics_pixels.md),
[`sc_metrics_clusters()`](https://jakubnowosad.com/supercells/reference/sc_metrics_clusters.md)

## Examples

``` r
library(supercells)
vol = terra::rast(system.file("raster/volcano.tif", package = "supercells"))
vol_sc = sc_slic(vol, step = 8, compactness = 1)
sc_metrics_global(vol, vol_sc)
#>   step compactness n_supercells mean_value_dist mean_spatial_dist
#> 1    8           1           90        2.007883          4.561378
#>   mean_combined_dist compactness_ratio_mean
#> 1           2.182844               3.532567
```
