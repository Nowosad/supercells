#' Pixel-level supercells metrics
#'
#' Computes per-pixel distance diagnostics from each pixel to its supercell center
#'
#' @param raster The input SpatRaster used to create `x`
#' @param x An sf object returned by [sc_slic()]
#' @param dist_fun A distance function name or function, as in [sc_slic()].
#' @param scale Logical. If `TRUE`, returns `spatial` and `value` as scaled
#' distances (`spatial_scaled`, `value_scaled`).
#' @param metrics Character vector of metrics to return. Options:
#' `"spatial"`, `"value"`, `"combined"`, `"balance"`.
#' Default: `c("spatial", "value", "combined")`.
#' @param compactness A compactness value used for the supercells
#' If missing, uses `attr(x, "compactness")` when available
#' @param step A step value used for the supercells
#' If missing, uses `attr(x, "step")` when available
#'
#' @details
#' If `x` lacks `supercells`, `x`, or `y` columns, they are derived from geometry
#' and row order, which may differ from the original centers.
#' When using SLIC0 (set `compactness = "auto"` in [sc_slic()]), combined and balance metrics use per-supercell
#' adaptive compactness (SLIC0), and scaled value distances are computed with the
#' per-supercell max value distance.
#'
#' @return A SpatRaster with one or more layers depending on `metrics`.
#' Interpretation:
#' \describe{
#'   \item{spatial}{Lower values indicate more compact supercells.}
#'   \item{value}{Lower values indicate more homogeneous supercells.}
#'   \item{combined}{Overall distance; mainly useful for ranking.}
#'   \item{balance}{0 indicates balance; negative values indicate spatial dominance;
#'   positive values indicate value dominance.}
#' }
#' Metrics:
#' \describe{
#'   \item{spatial}{Spatial distance from each pixel to its supercell center
#'   in grid-cell units (row/column index distance).}
#'   \item{value}{Value distance from each pixel to its supercell center in
#'   the raster value space.}
#'   \item{combined}{Combined distance using `compactness` and `step`.}
#'   \item{balance}{Signed log ratio of scaled value distance to scaled
#'   spatial distance; 0 indicates balance. Always computed from scaled components.}
#' }
#' When `scale = TRUE`, `spatial` and `value` are returned as
#' `spatial_scaled` and `value_scaled`.
#' @seealso [`sc_slic()`], [`sc_metrics_supercells()`], [`sc_metrics_global()`]
#' @export
#' @examples
#' library(supercells)
#' vol = terra::rast(system.file("raster/volcano.tif", package = "supercells"))
#' vol_sc = sc_slic(vol, step = 8, compactness = 7)
#' metrics = sc_metrics_pixels(vol, vol_sc, scale = TRUE)
#' terra::panel(metrics, nr = 1)
sc_metrics_pixels = function(raster, x, dist_fun = "euclidean", scale = TRUE,
                             metrics = c("spatial", "value", "combined"),
                             compactness, step) {
  
  if (any(!metrics %in% c("spatial", "value", "combined", "balance"))) {
    stop("metrics must be one or more of: spatial, value, combined, balance", call. = FALSE)
  }
  
  prep = .sc_metrics_prep(raster, x, dist_fun, compactness, step)

  out = sc_metrics_pixels_cpp(prep$clusters, prep$centers_xy, prep$centers_vals, prep$vals,
                              step = prep$step, compactness = prep$compactness,
                              adaptive_compactness = prep$adaptive_compactness,
                              dist_name = prep$dist_name, dist_fun = prep$dist_fun)

  spatial = terra::rast(out[["spatial"]])
  value = terra::rast(out[["value"]])
  combined = terra::rast(out[["combined"]])

  if (isTRUE(scale) || "balance" %in% metrics) {
    spatial_scaled = spatial / prep$step
    value_scaled = terra::rast(out[["value_scaled"]])
    if (isTRUE(scale)) {
      spatial = spatial_scaled
      value = value_scaled
    }
    if ("balance" %in% metrics) {
      balance = log(value_scaled / spatial_scaled)
    }
  } 

  if ("balance" %in% metrics) {
    result = c(spatial, value, combined, balance)
    names(result) = c("spatial", "value", "combined", "balance")
  } else {
    result = c(spatial, value, combined)
    names(result) = c("spatial", "value", "combined")
  }

  result = result[[metrics]]

  if (isTRUE(scale)) {
    names(result) = sub("^spatial$", "spatial_scaled", names(result))
    names(result) = sub("^value$", "value_scaled", names(result))
  }

  terra::ext(result) = terra::ext(prep$raster)
  terra::crs(result) = terra::crs(prep$raster)

  return(result)
}
