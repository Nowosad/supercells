#' Create supercells as points
#'
#' Runs the SLIC workflow and returns supercell centers as points. Use \code{iter = 0}
#' to return the initial centers (before iterations).
#' For polygon outputs, use `sc_slic`; for raster IDs, use `sc_slic_raster`.
#'
#' @inheritParams sc_slic
#' @return An sf object with supercell center points and summary statistics.
#' @export
sc_slic_points = function(x, k = NULL, step = NULL, centers = NULL, compactness,
                          dist_fun = "euclidean", avg_fun = "mean", clean = TRUE,
                          iter = 10, transform = NULL, minarea, metadata = FALSE,
                          chunks = FALSE, future = FALSE, verbose = 0,
                          iter_diagnostics = FALSE) {
  prep = .sc_slic_prep(x, k, step, centers, compactness, dist_fun, avg_fun,
                       minarea, chunks, iter, transform, metadata, iter_diagnostics)

  if (iter == 0) {
    clean = FALSE
  }

  points_sf = if (nrow(prep$chunk_ext) == 1) {
    list(.sc_slic_run_single_points(prep, compactness, clean, iter, transform, verbose, future))
  } else {
    .sc_slic_run_chunks_points(prep, compactness, clean, iter, transform, verbose, future)
  }

  iter_attr = NULL
  if (prep$iter_diagnostics && is.list(points_sf) && length(points_sf) > 0) {
    iter_attr = attr(points_sf[[1]], "iter_diagnostics")
  }
  .sc_slic_post(points_sf, metadata, prep$step, compactness, iter_attr)
}

.sc_slic_run_single_points = function(prep, compactness, clean, iter, transform, verbose, future) {
  ext = prep$chunk_ext[1, ]
  run_slic_chunk_points(ext, prep$x, step = prep$step, compactness = compactness,
                        dist_name = prep$funs$dist_name, dist_fun = prep$funs$dist_fun,
                        avg_fun_fun = prep$funs$avg_fun_fun, avg_fun_name = prep$funs$avg_fun_name,
                        clean = clean, iter = iter, minarea = prep$minarea, transform = transform,
                        input_centers = prep$input_centers, verbose = verbose,
                        iter_diagnostics = prep$iter_diagnostics)
}

.sc_slic_run_chunks_points = function(prep, compactness, clean, iter, transform, verbose, future) {
  .sc_slic_apply_chunks(prep, run_slic_chunk_points, compactness, clean, iter, transform, verbose, future)
}
