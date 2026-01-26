# supercells 1.7

* Added `sc_tune_compactness()` to estimate compactness from a light SLIC run and pixel metrics
* Metrics API update: `sc_metrics_clusters()` renamed to `sc_metrics_supercells()`; `metrics` now uses consistent names (`spatial`, `value`, `combined`, `balance`) across pixels/supercells/global, with scaled outputs suffixed by `_scaled`
* Balance metrics now use the absolute log ratio of scaled value/spatial distances for clearer interpretation
* `iter = 0` is now only supported for point outputs; `sc_slic()`/`sc_slic_raster()` error with guidance when used with polygons/rasters
* `transform` support was removed from `sc_slic()`, `sc_slic_points()`, and `sc_slic_raster()` to keep them minimal
* Legacy `supercells()` now uses `sc_slic()` while preserving `transform = "to_LAB"` behavior
* Verbose argument moved to the end of the argument list across R and C++ interfaces for consistency
* Empty-cluster handling is now consistent across connectivity cleanup, preserving prior centers/values

# supercells 1.6

* Spatial distance uses the more precise (fractional) center positions instead of rounding them down
* Neighborhood search uses rounded center positions, so results can differ slightly from 1.5
* Custom centers are handled more safely, with checks for invalid positions, etc.
* Local-minimum search is more robust at edges and avoids extra allocations
* Small-region cleanup is safer when no neighbor exists and avoids zero-division cases
* Empty-center cases are handled safely, keeping previous centers/values when needed
* Tests now use `tests/testthat/setup.R`, and expectations were updated for 1.6 behavior

# supercells 1.5

* Refactored core C++ SLIC implementation for clearer data flow, center updates, and diagnostics support, with related cleanup in supporting utilities
* Expanded internal R helper modules (`.sc_prep_raster()`, `.sc_metrics_prep()`, `.sc_slic_prep_args()`) to centralize validation, chunk handling, and function normalization across workflows
* Reorganized SLIC runners and chunk helpers (`run_slic_chunk_raster()`, `run_slic_chunk_points()`, `run_slic_chunks()`, `update_supercells_ids()`) to improve metadata handling, chunk ID normalization, and cross-output consistency
* Added and standardized SLIC functions `sc_slic()`, `sc_slic_points()`, and `sc_slic_raster()` with consistent options, metadata handling, and chunking behavior -- the `sc_slic()` is planned as the main user-facing function in the future (replacing `supercells()`)
* Added diagnostics API: `sc_metrics_pixels()`, `sc_metrics_clusters()`, and `sc_metrics_global()`, each documented with examples and output descriptions
* Added iteration diagnostics output via `iter_diagnostics` and a plotting helper `sc_plot_iter_diagnostics()` for quick inspection of convergence behavior
* Expanded tests for metrics and regression cases

# supercells 1.0.3

* C++ code improvements (Avoid implicit conversion from sexp to double, see #39)
* Improves documentation of the `fun_avg` argument (see #41)

# supercells 1.0.2

* Fixes centroid coordinates (see #32)

# supercells 1.0.0

* Added a `NEWS.md` file to track changes to the package.
* This version is a stable one as described in Nowosad and Stepinski (2022)
