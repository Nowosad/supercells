# Create supercells as a raster

Runs the SLIC workflow and returns a raster of supercell IDs IDs are
1-based and are unique across chunks when chunking is used For polygon
outputs, use
[`sc_slic()`](https://jakubnowosad.com/supercells/reference/sc_slic.md);
for point centers, use
[`sc_slic_points()`](https://jakubnowosad.com/supercells/reference/sc_slic_points.md)

## Usage

``` r
sc_slic_raster(
  x,
  step = NULL,
  compactness,
  dist_fun = "euclidean",
  avg_fun = "mean",
  clean = TRUE,
  minarea,
  iter = 10,
  k = NULL,
  centers = NULL,
  outcomes = "supercells",
  chunks = FALSE,
  verbose = 0
)
```

## Arguments

- x:

  An object of class SpatRaster (terra) or class stars (stars).

- step:

  Initial center spacing (alternative is `k`). Provide a plain numeric
  value for cell units, or use
  [`use_meters()`](https://jakubnowosad.com/supercells/reference/use_meters.md)
  for map-distance steps in meters (automatically converted to cells
  using raster resolution).

- compactness:

  A compactness value. Use
  [`sc_tune_compactness()`](https://jakubnowosad.com/supercells/reference/sc_tune_compactness.md)
  to estimate it. Use
  [`use_adaptive()`](https://jakubnowosad.com/supercells/reference/use_adaptive.md)
  to enable adaptive compactness (ASLIC).

- dist_fun:

  A supported distance name from
  [`sc_list_distances()`](https://jakubnowosad.com/supercells/reference/sc_list_distances.md)
  or a custom function that takes two numeric vectors and returns one
  finite, nonnegative value. Smaller values must mean greater
  similarity.

- avg_fun:

  An averaging function name or custom function used to summarize values
  within each supercell. Supported names: "mean" and "median". A custom
  function must accept a numeric vector and return a single numeric
  value.

- clean:

  Should connectivity of the supercells be enforced?

- minarea:

  When `clean = TRUE`, connected components containing `minarea` cells
  or fewer are merged into an adjacent component. If `minarea` is
  omitted or set to `0`, the threshold is selected automatically as one
  quarter of the raster cell count divided by the number of initialized
  centers: `floor(ncell / (4 * n_centers))`.

- iter:

  Number of iterations.

- k:

  The number of supercells desired (alternative to `step`).

- centers:

  Optional sf object of custom centers. Requires `step`.

- outcomes:

  Character vector controlling which fields are returned. Only
  `"supercells"` is supported in `sc_slic_raster()`.

- chunks:

  Chunking option. Use `FALSE` for no chunking, `TRUE` for automatic
  chunking based on size, or a numeric value for a fixed chunk size (in
  number of cells per side).

- verbose:

  Verbosity level.

## Value

A SpatRaster with supercell IDs.

## See also

[`sc_slic()`](https://jakubnowosad.com/supercells/reference/sc_slic.md),
[`sc_slic_points()`](https://jakubnowosad.com/supercells/reference/sc_slic_points.md)

## Examples

``` r
library(supercells)
vol = terra::rast(system.file("raster/volcano.tif", package = "supercells"))
vol_ids = sc_slic_raster(vol, step = 8, compactness = 1)
terra::plot(vol_ids)
```
