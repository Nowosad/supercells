# Cluster-level supercells metrics

Computes per-cluster distance diagnostics

## Usage

``` r
sc_metrics_clusters(raster, x, dist_fun = "euclidean", compactness, step)
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

An sf object with one row per supercell and columns:

- supercells:

  Supercell ID.

- mean_value_dist:

  Mean value distance from pixels to the supercell center in value
  space.

- mean_spatial_dist:

  Mean spatial distance from pixels to the supercell center in grid-cell
  units (row/column index distance).

- mean_combined_dist:

  Mean combined distance using `compactness` and `step` to scale value
  and spatial distances.

- compactness_ratio:

  Ratio of scaled value distance to scaled spatial distance; `NA` when
  `compactness` or `step` is zero.

## Details

If `x` lacks `supercells`, `x`, or `y` columns, they are derived from
geometry and row order, which may differ from the original centers

## See also

[`sc_slic()`](https://jakubnowosad.com/supercells/reference/sc_slic.md),
[`sc_metrics_pixels()`](https://jakubnowosad.com/supercells/reference/sc_metrics_pixels.md),
[`sc_metrics_global()`](https://jakubnowosad.com/supercells/reference/sc_metrics_global.md)

## Examples

``` r
library(supercells)
vol = terra::rast(system.file("raster/volcano.tif", package = "supercells"))
vol_sc = sc_slic(vol, step = 8, compactness = 1, metadata = TRUE)
cl = sc_metrics_clusters(vol, vol_sc)
head(cl)
#> Simple feature collection with 6 features and 5 fields
#> Geometry type: POLYGON
#> Dimension:     XY
#> Bounding box:  xmin: 2667400 ymin: 6479135 xmax: 2667550 ymax: 6479575
#> Projected CRS: NZGD49 / New Zealand Map Grid
#>   supercells mean_value_dist mean_spatial_dist mean_combined_dist
#> 1          1       0.2815283          3.658929          0.5642966
#> 2          2       0.2170883          4.446049          0.6287833
#> 3          3       0.7511050          4.540517          0.9926147
#> 4          4       0.4770408          3.313325          0.6410107
#> 5          5       0.8576279          4.953696          1.1359290
#> 6          6       0.9980867          5.395495          1.2809768
#>   compactness_ratio                       geometry
#> 1         0.6155426 POLYGON ((2667400 6479575, ...
#> 2         0.3906179 POLYGON ((2667470 6479515, ...
#> 3         1.3233825 POLYGON ((2667480 6479425, ...
#> 4         1.1518118 POLYGON ((2667450 6479375, ...
#> 5         1.3850312 POLYGON ((2667510 6479375, ...
#> 6         1.4798817 POLYGON ((2667500 6479305, ...
```
