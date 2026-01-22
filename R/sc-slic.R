#' Create supercells
#'
#' Creates supercells from single- or multi-band rasters using a SLIC-style algorithm.
#' The function supports either a target number of supercells (`k`) or a fixed grid
#' spacing (`step`), as well as optional custom centers and chunked processing.
#'
#' @details
#' Use `sc_slic` for polygon outputs. For raster IDs or point centers, see
#' `sc_slic_raster` and `sc_slic_points`.
#'
#' @param x An object of class SpatRaster (terra) or class stars (stars).
#' @param k The number of supercells desired (alternative to `step`).
#' @param step The distance (number of cells) between initial centers (alternative to `k`).
#' @param centers Optional sf object of custom centers. Requires `step`.
#' @param compactness A compactness value.
#' @param dist_fun A distance function name or a custom function. Supported names:
#' "euclidean", "jsd", "dtw", "dtw2d", or any method from `philentropy::getDistMethods()`.
#' A custom function must accept two numeric vectors and return a single numeric value.
#' @param avg_fun An averaging function name or custom function used to summarize
#' values within each supercell. Supported names: "mean" and "median". A custom
#' function must accept a numeric vector and return a single numeric value.
#' @param clean Should connectivity of the supercells be enforced?
#' @param iter Number of iterations.
#' @param transform Optional transformation applied before segmentation. Currently
#' supports "to_LAB" for RGB inputs.
#' @param minarea Minimal size of a supercell (in cells).
#' @param metadata Logical. Should metadata columns be kept?
#' @param chunks Chunking option. Use `FALSE` for no chunking, `TRUE` for
#' automatic chunking based on size, or a numeric value for a fixed chunk size
#' (in number of cells per side).
#' @param future Logical. Use future for parallelization?
#' @param verbose Verbosity level.
#' @param iter_diagnostics Logical. If `TRUE`, returns iteration diagnostics as an attribute
#' (`iter_diagnostics`) on the output. Only available when chunks are not used.
#'
#' @return An sf object with the supercell polygons and summary statistics.
#' Attributes `step` and `compactness` are attached to the result.
#' @export
sc_slic = function(x, k = NULL, step = NULL, centers = NULL, compactness, dist_fun = "euclidean",
                   avg_fun = "mean", clean = TRUE, iter = 10, transform = NULL,
                   minarea, metadata = FALSE, chunks = FALSE, future = FALSE, verbose = 0,
                   iter_diagnostics = FALSE) {
  prep = .sc_slic_prep(x, k, step, centers, compactness, dist_fun, avg_fun,
                       minarea, chunks, iter, transform, metadata, iter_diagnostics)

  if (iter == 0) {
    clean = FALSE
  }

  slic_sf = if (nrow(prep$chunk_ext) == 1) {
    .sc_slic_run_single_vector(prep, compactness, clean, iter, transform, verbose, future)
  } else {
    .sc_slic_run_chunks_vector(prep, compactness, clean, iter, transform, verbose, future)
  }

  iter_attr = NULL
  if (prep$iter_diagnostics && is.list(slic_sf) && length(slic_sf) > 0) {
    iter_attr = attr(slic_sf[[1]], "iter_diagnostics")
  }
  .sc_slic_post(slic_sf, metadata, prep$step, compactness, iter_attr)
}

.sc_slic_run_single_vector = function(prep, compactness, clean, iter, transform, verbose, future) {
  ext = prep$chunk_ext[1, ]
  run_slic_chunks(ext, prep$x, step = prep$step, compactness = compactness,
                  dist_name = prep$funs$dist_name, dist_fun = prep$funs$dist_fun,
                  avg_fun_fun = prep$funs$avg_fun_fun, avg_fun_name = prep$funs$avg_fun_name,
                  clean = clean, iter = iter, minarea = prep$minarea, transform = transform,
                  input_centers = prep$input_centers, verbose = verbose,
                  iter_diagnostics = prep$iter_diagnostics)
}

.sc_slic_run_chunks_vector = function(prep, compactness, clean, iter, transform, verbose, future) {
  .sc_slic_apply_chunks(prep, run_slic_chunks, compactness, clean, iter, transform, verbose, future)
}
