# SLIC convergence diagnostics

Runs SLIC and returns per-iteration mean combined distance. The output
can be plotted directly with
[`plot()`](https://rspatial.github.io/terra/reference/plot.html).

## Usage

``` r
sc_slic_convergence(
  x,
  step = NULL,
  compactness,
  dist_fun = "euclidean",
  avg_fun = "mean",
  clean = FALSE,
  minarea,
  iter = 10,
  k = NULL,
  centers = NULL,
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
  to enable adaptive compactness based on local feature-space distances.

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

  Connectivity cleaning is not applied to convergence diagnostics.
  Setting this to `TRUE` produces a warning and is otherwise ignored.

- minarea:

  Ignored because connectivity cleaning is not applied.

- iter:

  Number of iterations.

- k:

  The number of supercells desired (alternative to `step`).

- centers:

  Optional sf object of custom centers. Requires `step`.

- verbose:

  Verbosity level.

## Value

A data frame with class `sc_slic_convergence` and columns:

- iter:

  Iteration number.

- mean_distance:

  Mean combined distance across assigned cells at each iteration.

## See also

[`sc_slic()`](https://jakubnowosad.com/supercells/reference/sc_slic.md),
[`plot()`](https://rspatial.github.io/terra/reference/plot.html)

## Examples

``` r
library(supercells)
vol = terra::rast(system.file("raster/volcano.tif", package = "supercells"))
conv = sc_slic_convergence(vol, step = 8, compactness = 5, iter = 10)
plot(conv)
```
