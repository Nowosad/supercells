# Estimate a compactness value

Estimates a compactness value for a chosen raster scale. The current
implementation supports one tuning metric, `"local_variability"`, which
estimates compactness directly from local value variability.

## Usage

``` r
sc_tune_compactness(
  raster,
  step = NULL,
  dist_fun = "euclidean",
  metric = "local_variability",
  k = NULL,
  centers = NULL
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

- dist_fun:

  A distance function name or a custom function. Supported names:
  "euclidean", "jsd", "dtw", "dtw2d", or any method from
  [`philentropy::getDistMethods()`](https://drostlab.github.io/philentropy/reference/getDistMethods.html).
  A custom function must accept two numeric vectors and return a single
  numeric value.

- metric:

  Which compactness metric to return. Currently only
  `"local_variability"` is supported. The argument is kept for future
  additions. For `"local_variability"`,
  `compactness = median(local_mean(d_value)) / dim_scale`, where local
  means are computed in windows around initial centers and `dim_scale`
  is inferred from the number of raster layers and `dist_fun`. This
  keeps compactness adjustable to dimensionality without requiring a
  separate user-facing scaling argument.

- k:

  The number of supercells desired (alternative to `step`).

- centers:

  Optional sf object of custom initial centers. Requires `step`.

## Value

A one-row data frame with columns `step` and `compactness`.

## See also

[`sc_slic()`](https://jakubnowosad.com/supercells/reference/sc_slic.md),
[`use_meters()`](https://jakubnowosad.com/supercells/reference/use_meters.md),
[`use_adaptive()`](https://jakubnowosad.com/supercells/reference/use_adaptive.md)

## Examples

``` r
library(terra)
#> terra 1.9.1
vol = rast(system.file("raster/volcano.tif", package = "supercells"))
tune = sc_tune_compactness(vol, step = 8)
tune$compactness
#> [1] 8.990234
```
