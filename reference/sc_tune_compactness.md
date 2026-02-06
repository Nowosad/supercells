# Estimate compactness from a light SLIC run

Runs a short SLIC segmentation (default `iter = 1`) and uses cell-level
distances to estimate a compactness value where value and spatial
distances are balanced for the chosen `step`. The global estimate uses a
pixel-weighted median over the sampled cells, while the local estimate
uses a median of per-center mean distances.

## Usage

``` r
sc_tune_compactness(
  raster,
  step = NULL,
  step_unit = "cells",
  compactness = 1,
  metrics = "global",
  dist_fun = "euclidean",
  avg_fun = "mean",
  clean = TRUE,
  minarea,
  iter = 1,
  k = NULL,
  centers = NULL,
  sample_size = 10000,
  value_scale = "auto"
)
```

## Arguments

- raster:

  A `SpatRaster`.

- step:

  The distance (number of cells) between initial centers (alternative is
  `k`).

- step_unit:

  Units for `step`. Use "cells" for pixel units or "map" for map units.

- compactness:

  Starting compactness used for the initial short run.

- metrics:

  Which compactness metric to return: `"global"` or `"local"`. Default:
  `"global"`.

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

  Number of SLIC iterations for the pilot run.

- k:

  The number of supercells desired (alternative to `step`).

- centers:

  Optional sf object of custom centers. Requires `step`.

- sample_size:

  Optional limit on the number of pixels used for the summary (passed to
  [`terra::global()`](https://rspatial.github.io/terra/reference/global.html)
  as `maxcell`).

- value_scale:

  Optional scale factor applied to the median value distance before
  computing compactness. Use `"auto"` to divide by `sqrt(nlyr(raster))`
  (useful for high-dimensional embeddings). Default: `"auto"`.

## Value

A one-row data frame with columns `step`, `metric`, and `compactness`.

## See also

[`sc_slic()`](https://jakubnowosad.com/supercells/reference/sc_slic.md)

## Examples

``` r
library(terra)
#> terra 1.8.93
vol = rast(system.file("raster/volcano.tif", package = "supercells"))
tune = sc_tune_compactness(vol, step = 8)
tune$compactness
#> [1] 6.864497
```
