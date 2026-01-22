#' Create supercells as points
#'
#' Runs the SLIC workflow and returns supercell centroids as points.
#' Use \code{iter = 0} to return the initial centers before iterations.
#' For polygon outputs, use `sc_slic`; for raster IDs, use `sc_slic_raster`
#'
#' @inheritParams sc_slic
#' @seealso [`sc_slic()`]
#' @return An sf object with supercell centroid points and summary statistics
#'
#' @export
#'
#' @examples
#' library(supercells)
#' vol = terra::rast(system.file("raster/volcano.tif", package = "supercells"))
#'
#' init_pts = sc_slic_points(vol, step = 20, compactness = 1, iter = 1)
#' terra::plot(vol)
#' plot(sf::st_geometry(init_pts), add = TRUE, pch = 3, col = "red")
#'
#' vol_pts = sc_slic_points(vol, step = 8, compactness = 1)
#' terra::plot(vol)
#' plot(sf::st_geometry(vol_pts), add = TRUE, pch = 16, cex = 0.4)
sc_slic_points = function(x, step = NULL, compactness, dist_fun = "euclidean",
                          avg_fun = "mean", clean = TRUE, minarea, iter = 10,
                          transform = NULL, k = NULL, centers = NULL,
                          metadata = FALSE, chunks = FALSE, future = FALSE, verbose = 0,
                          iter_diagnostics = FALSE) {
  prep_args = .sc_slic_prep_args(x, step, compactness, k, centers, dist_fun, avg_fun,
                            minarea, chunks, iter, transform, metadata, iter_diagnostics)

  if (iter == 0) {
    clean = FALSE
  }

  points_sf = if (nrow(prep_args$chunk_ext) == 1) {
    .sc_slic_run_single_points(prep_args, compactness, clean, iter, transform, verbose, future, metadata)
  } else {
    .sc_slic_run_chunks_points(prep_args, compactness, clean, iter, transform, verbose, future, metadata)
  }

  iter_attr = NULL
  if (prep_args$iter_diagnostics && is.list(points_sf) && length(points_sf) > 0) {
    iter_attr = attr(points_sf[[1]], "iter_diagnostics")
  }

  results = .sc_slic_post(points_sf, metadata, prep_args$step, compactness, iter_attr)
  return(results)
}

.sc_slic_run_single_points = function(prep, compactness, clean, iter, transform, verbose, future, metadata) {
  ext = prep$chunk_ext[1, ]
  results = list(run_slic_chunk_points(ext, prep$x, step = prep$step, compactness = compactness,
                             dist_name = prep$funs$dist_name, dist_fun = prep$funs$dist_fun,
                             avg_fun_fun = prep$funs$avg_fun_fun, avg_fun_name = prep$funs$avg_fun_name,
                             clean = clean, iter = iter, minarea = prep$minarea, transform = transform,
                             input_centers = prep$input_centers, verbose = verbose,
                             iter_diagnostics = prep$iter_diagnostics, metadata = metadata))
  return(results)
}

.sc_slic_run_chunks_points = function(prep, compactness, clean, iter, transform, verbose, future, metadata) {
  results = .sc_slic_apply_chunks(prep, run_slic_chunk_points, compactness, clean, iter, transform, verbose, future,
                                  metadata = metadata)
  return(results)
}
