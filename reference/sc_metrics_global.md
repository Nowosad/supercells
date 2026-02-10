# Global supercells metrics

Computes global distance diagnostics for supercells

## Usage

``` r
sc_metrics_global(
  x,
  sc,
  metrics = c("spatial", "value", "combined", "balance"),
  scale = TRUE,
  step,
  compactness,
  dist_fun = NULL
)
```

## Arguments

- x:

  The input SpatRaster used to create `sc`.

- sc:

  An sf object returned by
  [`sc_slic()`](https://jakubnowosad.com/supercells/reference/sc_slic.md).

- metrics:

  Character vector of metric ideas to return. Options: `"spatial"`,
  `"value"`, `"combined"`, `"balance"`. Default:
  `c("spatial", "value", "combined", "balance")`.

- scale:

  Logical. If `TRUE`, scales spatial and value distances; output columns
  are named with the `_scaled` suffix.

- step:

  A step value used for the supercells If missing, uses
  `attr(sc, "step")` when available

- compactness:

  A compactness value used for the supercells If missing, uses
  `attr(sc, "compactness")` when available. Compactness mode is read
  from `attr(sc, "compactness_method")` when available.

- dist_fun:

  A distance function name or function, as in
  [`sc_slic()`](https://jakubnowosad.com/supercells/reference/sc_slic.md).
  If missing or `NULL`, uses `attr(sc, "dist_fun")` when available.

## Value

A data.frame with a single row and columns:

- step:

  Step size used to generate supercells. Returned in meters when the
  input used `step = use_meters(...)`, otherwise in cells.

- compactness:

  Compactness value used to generate supercells; `NA` for adaptive
  compactness.

- compactness_method:

  Compactness method: `"constant"` for fixed compactness, `"local_max"`
  for adaptive compactness.

- n_supercells:

  Number of supercells with at least one non-missing pixel.

- mean_value_dist / mean_value_dist_scaled:

  Mean per-supercell value distance from cells to their supercell
  centers, averaged across supercells. Returned as `mean_value_dist` (or
  `mean_value_dist_scaled` when `scale = TRUE`). Lower values indicate
  more homogeneous supercells.

- mean_spatial_dist / mean_spatial_dist_scaled:

  Mean per-supercell spatial distance from cells to their supercell
  centers, averaged across supercells; units are grid cells (row/column
  index distance). If the input supercells were created with
  `step = use_meters(...)`, distances are reported in meters. Returned
  as `mean_spatial_dist` (or `mean_spatial_dist_scaled` when
  `scale = TRUE`). Lower values indicate more compact supercells.

- mean_combined_dist:

  Mean per-supercell combined distance, computed from value and spatial
  distances using `compactness` and `step`, averaged across supercells.
  Returned as `mean_combined_dist`. Lower values indicate lower overall
  distance and are mainly useful for ranking.

- balance:

  Mean signed log ratio of scaled value distance to scaled spatial
  distance (averaged across supercells); 0 indicates balance between
  value and spatial terms, negative values indicate spatial dominance,
  and positive values indicate value dominance.

## Details

Requires `sc` with metadata columns (`supercells`, `x`, `y`). If they
are missing, they are derived from geometry and row order. Use
`outcomes = c("supercells", "coordinates", "values")` when calling
[`sc_slic()`](https://jakubnowosad.com/supercells/reference/sc_slic.md)
or
[`supercells()`](https://jakubnowosad.com/supercells/reference/supercells.md)
to preserve original centers and IDs. Metrics are averaged across
supercells (each supercell has equal weight). When using SLIC0 (set
`compactness = use_adaptive()` in
[`sc_slic()`](https://jakubnowosad.com/supercells/reference/sc_slic.md)),
combined and balance metrics use per-supercell adaptive compactness
(SLIC0), and scaled value distances are computed with the per-supercell
max value distance.

## See also

[`sc_slic()`](https://jakubnowosad.com/supercells/reference/sc_slic.md),
[`sc_metrics_pixels()`](https://jakubnowosad.com/supercells/reference/sc_metrics_pixels.md),
[`sc_metrics_supercells()`](https://jakubnowosad.com/supercells/reference/sc_metrics_supercells.md)

## Examples

``` r
library(supercells)
vol = terra::rast(system.file("raster/volcano.tif", package = "supercells"))
vol_sc = sc_slic(vol, step = 8, compactness = 7)
sc_metrics_global(vol, vol_sc)
#>   step compactness compactness_method n_supercells mean_spatial_dist_scaled
#> 1    8           7           constant           88                0.4718607
#>   mean_value_dist_scaled mean_combined_dist    balance
#> 1              0.3701397          0.6517259 -0.3367309
```
