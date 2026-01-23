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
