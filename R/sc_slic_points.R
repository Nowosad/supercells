#' Create supercells as points
#'
#' Runs the SLIC workflow and returns supercell centers as points.
#' Use \code{iter = 0} to return the initial centers before iterations.
#' For polygon outputs, use [`sc_slic()`]; for raster output, use [`sc_slic_raster()`]
#' By default, only value summaries are returned; add
#' `outcomes = c("supercells", "coordinates", "values")` to include ids and x/y.
#'
#' @inheritParams sc_slic
#' @seealso [`sc_slic()`], [`sc_slic_raster()`]
#' @return An sf object with supercell center points and summary statistics
#'
#' @export
#'
#' @examples
#' library(supercells)
#' vol = terra::rast(system.file("raster/volcano.tif", package = "supercells"))
#' 
#' # initial centers (only after local minima placement, no iterations)
#' init_pts = sc_slic_points(vol, step = 12, compactness = 1, iter = 0)
#' terra::plot(vol)
#' plot(sf::st_geometry(init_pts), add = TRUE, pch = 3, col = "red")
#'
#' # final supercell centers
#' vol_pts = sc_slic_points(vol, step = 12, compactness = 1)
#' terra::plot(vol)
#' plot(sf::st_geometry(vol_pts), add = TRUE, pch = 16, col = "red")
sc_slic_points = function(x, step = NULL, compactness, dist_fun = "euclidean",
                          avg_fun = "mean", clean = TRUE, minarea, iter = 10,
                          k = NULL, centers = NULL,
                          outcomes = "values", chunks = FALSE,
                          verbose = 0) {
  if (iter == 0) {
    clean = FALSE
  }
  prep_args = .sc_slic_prep_args(x, step, compactness, dist_fun, avg_fun, clean, minarea, iter,
                                k, centers, outcomes, chunks, verbose)

  segment = .sc_slic_segment(prep_args, .sc_run_full_points, .sc_run_chunk_points)

  results = .sc_slic_post(segment$chunks, prep_args)
  return(results)
}
