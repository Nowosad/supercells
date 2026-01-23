#' Pixel-level supercells metrics
#'
#' Computes per-pixel spatial, value, and combined distance diagnostics
#'
#' @param raster The input SpatRaster used to create `x`
#' @param x An sf object returned by [sc_slic()]
#' @param dist_fun A distance function name or function, as in [sc_slic()]
#' @param compactness A compactness value used for the supercells
#' If missing, uses `attr(x, "compactness")` when available
#' @param step A step value used for the supercells
#' If missing, uses `attr(x, "step")` when available
#'
#' @details
#' If `x` lacks `supercells`, `x`, or `y` columns, they are derived from geometry
#' and row order, which may differ from the original centers
#'
#' @return A SpatRaster with three layers:
#' \describe{
#'   \item{spatial}{Spatial distance from each pixel to its supercell center
#'   in grid-cell units (row/column index distance).}
#'   \item{value}{Value distance from each pixel to its supercell center in
#'   the raster value space.}
#'   \item{combined}{Combined distance using `compactness` and `step` to
#'   scale value and spatial distances.}
#' }
#' @seealso [`sc_slic()`], [`sc_metrics_clusters()`], [`sc_metrics_global()`]
#' @export
#' @examples
#' library(supercells)
#' vol = terra::rast(system.file("raster/volcano.tif", package = "supercells"))
#' vol_sc = sc_slic(vol, step = 8, compactness = 1)
#' metrics = sc_metrics_pixels(vol, vol_sc)
#' terra::panel(metrics, nr = 1)
sc_metrics_pixels = function(raster, x, dist_fun = "euclidean", compactness, step) {
  prep = .sc_metrics_prep(raster, x, dist_fun, compactness, step)

  out = sc_metrics_pixels_cpp(prep$clusters, prep$centers_xy, prep$centers_vals, prep$vals,
                              step = prep$step, compactness = prep$compactness,
                              dist_name = prep$dist_name, dist_fun = prep$dist_fun)

  spatial = terra::rast(out[["spatial"]])
  value = terra::rast(out[["value"]])
  combined = terra::rast(out[["combined"]])

  result = c(spatial, value, combined)
  terra::ext(result) = terra::ext(prep$raster)
  terra::crs(result) = terra::crs(prep$raster)
  names(result) = c("spatial", "value", "combined")
  return(result)
}
