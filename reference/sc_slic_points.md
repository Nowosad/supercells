# Create supercells as points

Runs the SLIC workflow and returns supercell centers as points. Use
`iter = 0` to return the initial centers before iterations. For polygon
outputs, use
[`sc_slic()`](https://jakubnowosad.com/supercells/reference/sc_slic.md);
for raster output, use
[`sc_slic_raster()`](https://jakubnowosad.com/supercells/reference/sc_slic_raster.md)
By default, only value summaries are returned; add
`outcomes = c("supercells", "coordinates", "values")` to include ids and
x/y.

## Usage

``` r
sc_slic_points(
  x,
  step = NULL,
  compactness,
  dist_fun = "euclidean",
  avg_fun = "mean",
  clean = TRUE,
  minarea,
  iter = 10,
  step_unit = "cells",
  k = NULL,
  centers = NULL,
  outcomes = "values",
  chunks = FALSE,
  iter_diagnostics = FALSE,
  verbose = 0
)
```

## Arguments

- x:

  An object of class SpatRaster (terra) or class stars (stars).

- step:

  The distance (number of cells) between initial centers (alternative is
  `k`).

- compactness:

  A compactness value. Use
  [`sc_tune_compactness()`](https://jakubnowosad.com/supercells/reference/sc_tune_compactness.md)
  to estimate it. Set `compactness = "auto"` to enable adaptive
  compactness (SLIC0).

- dist_fun:

  A distance function name or a custom function. Supported names:
  "euclidean", "jsd", "dtw", "dtw2d", or any method from
  [`philentropy::getDistMethods()`](https://drostlab.github.io/philentropy/reference/getDistMethods.html).
  A custom function must accept two numeric vectors and return a single
  numeric value.

- avg_fun:

  An averaging function name or custom function used to summarize values
  within each supercell. Supported names: "mean" and "median". A custom
  function must accept a numeric vector and return a single numeric
  value.

- clean:

  Should connectivity of the supercells be enforced?

- minarea:

  Minimal size of a supercell (in cells).

- iter:

  Number of iterations.

- step_unit:

  Units for `step`. Use "cells" for pixel units or "map" for map units
  (converted to cells using raster resolution).

- k:

  The number of supercells desired (alternative to `step`).

- centers:

  Optional sf object of custom centers. Requires `step`.

- outcomes:

  Character vector controlling which fields are returned. Allowed values
  are "supercells", "coordinates", and "values". Default is "values".
  Use `outcomes = c("supercells", "coordinates", "values")` for full
  output.

- chunks:

  Chunking option. Use `FALSE` for no chunking, `TRUE` for automatic
  chunking based on size, or a numeric value for a fixed chunk size (in
  number of cells per side).

- iter_diagnostics:

  Logical. If `TRUE`, attaches iteration diagnostics as an attribute
  (`iter_diagnostics`) on the output. Only available when chunks are not
  used.

- verbose:

  Verbosity level.

## Value

An sf object with supercell center points and summary statistics

## See also

[`sc_slic()`](https://jakubnowosad.com/supercells/reference/sc_slic.md),
[`sc_slic_raster()`](https://jakubnowosad.com/supercells/reference/sc_slic_raster.md)

## Examples

``` r
library(supercells)
vol = terra::rast(system.file("raster/volcano.tif", package = "supercells"))

# initial centers (only after local minima placement, no iterations)
init_pts = sc_slic_points(vol, step = 12, compactness = 1, iter = 0)
terra::plot(vol)
plot(sf::st_geometry(init_pts), add = TRUE, pch = 3, col = "red")


# final supercell centers
vol_pts = sc_slic_points(vol, step = 12, compactness = 1)
terra::plot(vol)
plot(sf::st_geometry(vol_pts), add = TRUE, pch = 16, col = "red")
```
