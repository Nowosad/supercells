# Create supercells

Creates supercells from single- or multi-band rasters using an extended
SLIC algorithm. The function supports either a target number of
supercells (`k`) or a fixed grid spacing (`step`), as well as optional
custom centers and chunked processing.

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
  transform = NULL,
  k = NULL,
  centers = NULL,
  metadata = FALSE,
  chunks = FALSE,
  future = FALSE,
  verbose = 0,
  iter_diagnostics = FALSE
)
```

## Arguments

- x:

  An object of class SpatRaster (terra) or class stars (stars).

- step:

  The distance (number of cells) between initial centers (alternative to
  `k`).

- compactness:

  A compactness value.

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

- transform:

  Optional transformation applied before segmentation. Currently
  supports "to_LAB" for RGB inputs.

- k:

  The number of supercells desired (alternative to `step`).

- centers:

  Optional sf object of custom centers. Requires `step`.

- metadata:

  Logical. Should metadata columns be kept?

- chunks:

  Chunking option. Use `FALSE` for no chunking, `TRUE` for automatic
  chunking based on size, or a numeric value for a fixed chunk size (in
  number of cells per side).

- future:

  Logical. Use future for parallelization?

- verbose:

  Verbosity level.

- iter_diagnostics:

  Logical. If `TRUE`, returns iteration diagnostics as an attribute
  (`iter_diagnostics`) on the output. Only available when chunks are not
  used.

## Value

An sf object with the supercell polygons and summary statistics.
Information on `step` and `compactness` are attached to the result as
attributes.

## Details

Use `sc_slic` for polygon outputs. For raster or point centers outputs,
see `sc_slic_raster` and `sc_slic_points`. Evaluation and diagnostic
options:

- Iteration diagnostics: set `iter_diagnostics = TRUE` to attach an
  `iter_diagnostics` attribute (only available without chunking).

- Pixel diagnostics:
  [`sc_metrics_pixels()`](https://jakubnowosad.com/supercells/reference/sc_metrics_pixels.md)
  for per-pixel spatial, value, and combined distances.

- Cluster diagnostics:
  [`sc_metrics_clusters()`](https://jakubnowosad.com/supercells/reference/sc_metrics_clusters.md)
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

[`sc_slic_raster()`](https://jakubnowosad.com/supercells/reference/sc_slic_raster.md),
[`sc_slic_points()`](https://jakubnowosad.com/supercells/reference/sc_slic_points.md),
[`sc_plot_iter_diagnostics()`](https://jakubnowosad.com/supercells/reference/sc_plot_iter_diagnostics.md),
[`sc_metrics_pixels()`](https://jakubnowosad.com/supercells/reference/sc_metrics_pixels.md),
[`sc_metrics_clusters()`](https://jakubnowosad.com/supercells/reference/sc_metrics_clusters.md),
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
