# supercells 1.9

* Added `outcomes` argument to `sc_slic()`, `sc_slic_points()`, and `sc_slic_raster()`; replaces `metadata` for controlling returned fields
* Added experimental `sc_merge_supercells()` for adjacency-constrained greedy merging
* Added `sc_dist_vec_cpp()` (C++ distance wrapper) to support merge utilities
* Documentation and vignettes updated (pkgdown refresh, new articles, and revised examples)

# supercells 1.8

* `compactness = "auto"` enables SLIC0-style adaptive compactness
* Metrics now support SLIC0 adaptive compactness (pixels/supercells/global)
* `sc_tune_compactness()` returns only `step` + one compactness, selected by `metrics = "global"` or `"local"`
* Chunked raster IDs now use a consistent center-count offset strategy

# supercells 1.7

* Added `sc_tune_compactness()` to estimate compactness from a short SLIC run
* Metrics API cleanup: renamed `sc_metrics_clusters()` to `sc_metrics_supercells()`, unified metric names, `_scaled` suffix for scaled outputs
* Balance metrics now use absolute log ratio of scaled value/spatial distances
* `iter = 0` only works for point outputs; polygons/rasters error with guidance
* `sc_slic()`/`sc_slic_points()`/`sc_slic_raster()` dropped the `transform` argument; legacy `supercells()` keeps `transform = "to_LAB"`
* New `step_unit` supports map-unit step sizes
* Chunking improvements: conservative memory estimates, `options(supercells.chunk_mem_gb)`, chunk sizes align to `step`, deterministic ID offsets with file-backed merge
* Verbose argument moved to the end in R and C++ APIs
* Empty-cluster handling now consistent across cleanup
* Removed `future`-based parallel chunking

# supercells 1.6

* Spatial distance uses precise (fractional) center positions; neighborhood search still uses rounded centers (minor output differences are expected vs 1.5)
* More robust center handling (custom centers, empty centers) and safer cleanup for small regions
* Local-minimum search improved
* Tests standardized via `tests/testthat/setup.R` and updated expectations

# supercells 1.5

* Major C++ SLIC refactor (clearer data flow, diagnostics support) with supporting utility cleanup
* Standardized `sc_slic()`, `sc_slic_points()`, `sc_slic_raster()` with consistent options/metadata/chunking; `sc_slic()` is now the main entry point
* Added diagnostics API (`sc_metrics_pixels()`, `sc_metrics_clusters()`, `sc_metrics_global()`) and `iter_diagnostics` + `sc_plot_iter_diagnostics()`
* Centralized internal helpers for validation, chunking, and normalization
* Expanded regression and metrics tests

# supercells 1.0.3

* C++ code improvements (Avoid implicit conversion from sexp to double, see #39)
* Improves documentation of the `fun_avg` argument (see #41)

# supercells 1.0.2

* Fixes centroid coordinates (see #32)

# supercells 1.0.0

* Added a `NEWS.md` file to track changes to the package.
* This version is a stable one as described in Nowosad and Stepinski (2022)
