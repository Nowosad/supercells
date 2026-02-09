# Create supercells

Creates supercells from single- or multi-band rasters using an extended
SLIC algorithm. The function supports either a target number of
supercells (`k`) or a fixed grid spacing (`step`), as well as optional
custom centers and chunking.

## Usage

``` r
sc_slic(
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
  outcomes = c("supercells", "coordinates", "values"),
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
  to enable adaptive compactness (SLIC0).

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

- k:

  The number of supercells desired (alternative to `step`).

- centers:

  Optional sf object of custom centers. Requires `step`.

- outcomes:

  Character vector controlling which fields are returned. Allowed values
  are "supercells", "coordinates", and "values". Default is full output
  (`c("supercells", "coordinates", "values")`). Use
  `outcomes = "values"` for value summaries only.

- chunks:

  Chunking option. Use `FALSE` for no chunking, `TRUE` for automatic
  chunking based on size, or a numeric value for a fixed chunk size (in
  number of cells per side).

- verbose:

  Verbosity level.

## Value

An sf object with the supercell polygons and summary statistics.
Information on `step`, `compactness`, and `adaptive_method` are attached
to the result as attributes (`adaptive_method` is `NULL` for fixed
compactness).

## Details

Use `sc_slic()` for polygon outputs. For raster or point centers
outputs, see
[`sc_slic_raster()`](https://jakubnowosad.com/supercells/reference/sc_slic_raster.md)
and
[`sc_slic_points()`](https://jakubnowosad.com/supercells/reference/sc_slic_points.md).
Evaluation and diagnostic options:

- Iteration convergence: use
  [`sc_slic_convergence()`](https://jakubnowosad.com/supercells/reference/sc_slic_convergence.md)
  and plot its output.

- Pixel diagnostics:
  [`sc_metrics_pixels()`](https://jakubnowosad.com/supercells/reference/sc_metrics_pixels.md)
  for per-pixel spatial, value, and combined distances.

- Cluster diagnostics:
  [`sc_metrics_supercells()`](https://jakubnowosad.com/supercells/reference/sc_metrics_supercells.md)
  for per-supercell summaries.

- Global diagnostics:
  [`sc_metrics_global()`](https://jakubnowosad.com/supercells/reference/sc_metrics_global.md)
  for a single-row summary.

## References

Achanta, R., Shaji, A., Smith, K., Lucchi, A., Fua, P., & Süsstrunk, S.
(2012). SLIC Superpixels Compared to State-of-the-Art Superpixel
Methods. IEEE Transactions on Pattern Analysis and Machine Intelligence,
34(11), 2274–2282. https://doi.org/10.1109/tpami.2012.120

Nowosad, J., Stepinski, T. (2022). Extended SLIC superpixels algorithm
for applications to non-imagery geospatial rasters. International
Journal of Applied Earth Observation and Geoinformation,
https://doi.org/10.1016/j.jag.2022.102935

## See also

[`use_meters()`](https://jakubnowosad.com/supercells/reference/use_meters.md),
[`use_adaptive()`](https://jakubnowosad.com/supercells/reference/use_adaptive.md),
[`sc_slic_raster()`](https://jakubnowosad.com/supercells/reference/sc_slic_raster.md),
[`sc_slic_points()`](https://jakubnowosad.com/supercells/reference/sc_slic_points.md),
[`sc_slic_convergence()`](https://jakubnowosad.com/supercells/reference/sc_slic_convergence.md),
[`sc_metrics_pixels()`](https://jakubnowosad.com/supercells/reference/sc_metrics_pixels.md),
[`sc_metrics_supercells()`](https://jakubnowosad.com/supercells/reference/sc_metrics_supercells.md),
[`sc_metrics_global()`](https://jakubnowosad.com/supercells/reference/sc_metrics_global.md)

## Examples

``` r
library(supercells)
# One variable

vol = terra::rast(system.file("raster/volcano.tif", package = "supercells"))
vol_slic1 = sc_slic(vol, step = 8, compactness = 1)
terra::plot(vol)
plot(sf::st_geometry(vol_slic1), add = TRUE, lwd = 0.2)
```
