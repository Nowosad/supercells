#' Global supercells metrics
#'
#' Computes global distance diagnostics for supercells
#'
#' @details
#' Requires `x` with metadata columns (`supercells`, `x`, `y`)
#' If they are missing, they are derived from geometry and row order
#' Set `metadata = TRUE` when calling `sc_slic()` or `supercells()`
#'
#' @inheritParams sc_metrics_pixels
#' @return A data.frame with a single row of global metrics
#' @export
#' @examples
#' library(supercells)
#' vol = terra::rast(system.file("raster/volcano.tif", package = "supercells"))
#' vol_sc = sc_slic(vol, step = 8, compactness = 1)
#' sc_metrics_global(vol, vol_sc)
sc_metrics_global = function(raster, x, dist_fun = "euclidean", compactness, step) {
  prep = .sc_metrics_prep(raster, x, dist_fun, compactness, step)
  out = sc_metrics_global_cpp(prep$clusters, prep$centers_xy, prep$centers_vals, prep$vals,
                              step = prep$step, compactness = prep$compactness,
                              dist_name = prep$dist_name, dist_fun = prep$dist_fun)
  
  results = cbind(data.frame(step = prep$step, compactness = prep$compactness), out)
  return(results)
}
