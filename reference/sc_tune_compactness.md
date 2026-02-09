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
  compactness = 1,
  metric = "global",
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

  Initial center spacing (alternative is `k`). Provide a plain numeric
  value for cell units, or use
  [`use_meters()`](https://jakubnowosad.com/supercells/reference/use_meters.md)
  for map-distance steps in meters (automatically converted to cells
  using raster resolution).

- compactness:

  Starting compactness used for the initial short run.

- metric:

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

  Scale factor for value distances during tuning. Global metric:
  `compactness = (median(d_value) / value_scale) * step / median(d_spatial)`.
  Local metric:
  `compactness = median(local_mean(d_value) / value_scale)`. `"auto"`
  uses `sqrt(nlyr(raster))` (good for Euclidean-like distances); for
  bounded/angular distances (e.g., cosine), `value_scale = 1` is often
  better. Default: `"auto"`.

## Value

A one-row data frame with columns `step`, `metric`, `dist_fun`, and
`compactness`.

## See also

[`sc_slic()`](https://jakubnowosad.com/supercells/reference/sc_slic.md),
[`use_meters()`](https://jakubnowosad.com/supercells/reference/use_meters.md),
[`use_adaptive()`](https://jakubnowosad.com/supercells/reference/use_adaptive.md)

## Examples

``` r
library(terra)
#> terra 1.8.93
vol = rast(system.file("raster/volcano.tif", package = "supercells"))
tune = sc_tune_compactness(vol, step = 8)
tune$compactness
#> [1] 6.864497
```
