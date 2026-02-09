# Supercell-level metrics

Computes per-supercell distance diagnostics

## Usage

``` r
sc_metrics_supercells(
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
  `attr(sc, "compactness")` when available. Adaptive mode is read from
  `attr(sc, "adaptive_method")` when available.

- dist_fun:

  A distance function name or function, as in
  [`sc_slic()`](https://jakubnowosad.com/supercells/reference/sc_slic.md).
  If missing or `NULL`, uses `attr(sc, "dist_fun")` when available.

## Value

An sf object with one row per supercell and columns: Interpretation:

- mean_value_dist:

  Lower values indicate more homogeneous supercells.

- mean_spatial_dist:

  Lower values indicate more compact supercells.

- mean_combined_dist:

  Overall distance; mainly useful for ranking.

- balance:

  0 indicates balance between value and spatial terms; negative values
  indicate spatial dominance; positive values indicate value dominance.

Metrics:

- supercells:

  Supercell ID.

- spatial:

  Mean spatial distance from cells to the supercell center in grid-cell
  units (row/column index distance). If the input supercells were
  created with `step = use_meters(...)`, distances are reported in
  meters. Returned as `mean_spatial_dist` (or `mean_spatial_dist_scaled`
  when `scale = TRUE`).

- value:

  Mean value distance from cells to the supercell center in value space.
  Returned as `mean_value_dist` (or `mean_value_dist_scaled` when
  `scale = TRUE`).

- combined:

  Mean combined distance using `compactness` and `step` to scale value
  and spatial distances. Returned as `mean_combined_dist`.

- balance:

  Signed log ratio of scaled value distance to scaled spatial distance;
  0 indicates balance.

## Details

If `sc` lacks `supercells`, `x`, or `y` columns, they are derived from
geometry and row order, which may differ from the original centers When
using SLIC0 (set `compactness = use_adaptive()` in
[`sc_slic()`](https://jakubnowosad.com/supercells/reference/sc_slic.md)),
combined and balance metrics use per-supercell adaptive compactness
(SLIC0), and scaled value distances are computed with the per-supercell
max value distance.

## See also

[`sc_slic()`](https://jakubnowosad.com/supercells/reference/sc_slic.md),
[`sc_metrics_pixels()`](https://jakubnowosad.com/supercells/reference/sc_metrics_pixels.md),
[`sc_metrics_global()`](https://jakubnowosad.com/supercells/reference/sc_metrics_global.md)

## Examples

``` r
library(supercells)
vol = terra::rast(system.file("raster/volcano.tif", package = "supercells"))
vol_sc = sc_slic(vol, step = 8, compactness = 7)
cl = sc_metrics_supercells(vol, vol_sc)
head(cl)
#> Simple feature collection with 6 features and 5 fields
#> Geometry type: POLYGON
#> Dimension:     XY
#> Bounding box:  xmin: 2667400 ymin: 6479045 xmax: 2667500 ymax: 6479575
#> Projected CRS: NZGD49 / New Zealand Map Grid
#>   supercells mean_spatial_dist_scaled mean_value_dist_scaled mean_combined_dist
#> 1          1                0.4195350             0.08554244          0.4325501
#> 2          2                0.4539066             0.04427186          0.4587638
#> 3          3                0.4003193             0.14495014          0.4452615
#> 4          4                0.4570244             0.14389321          0.4939813
#> 5          5                0.4450208             0.31268263          0.5773152
#> 6          6                0.4680141             0.20287698          0.5398852
#>      balance                       geometry
#> 1 -1.5901344 POLYGON ((2667400 6479575, ...
#> 2 -2.3275422 POLYGON ((2667400 6479495, ...
#> 3 -1.0158726 POLYGON ((2667440 6479415, ...
#> 4 -1.1556654 POLYGON ((2667460 6479345, ...
#> 5 -0.3529323 POLYGON ((2667460 6479265, ...
#> 6 -0.8358987 POLYGON ((2667450 6479175, ...
```
