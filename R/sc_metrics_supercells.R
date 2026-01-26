#' Supercell-level metrics
#'
#' Computes per-supercell distance diagnostics
#'
#' @inheritParams sc_metrics_pixels
#' @param metrics Character vector of metric ideas to return. Options:
#' `"spatial"`, `"value"`, `"combined"`, `"balance"`. Default:
#' `c("spatial", "value", "combined", "balance")`.
#' @details
#' If `x` lacks `supercells`, `x`, or `y` columns, they are derived from geometry
#' and row order, which may differ from the original centers
#' @return An sf object with one row per supercell and columns:
#' Interpretation:
#' \describe{
#'   \item{mean_value_dist}{Lower values indicate more homogeneous supercells.}
#'   \item{mean_spatial_dist}{Lower values indicate more compact supercells.}
#'   \item{mean_combined_dist}{Overall distance; mainly useful for ranking.}
#'   \item{balance}{0 indicates balance between value and spatial terms.}
#' }
#' Metrics:
#' \describe{
#'   \item{supercells}{Supercell ID.}
#'   \item{spatial}{Mean spatial distance from pixels to the supercell center
#'   in grid-cell units (row/column index distance). Returned as
#'   `mean_spatial_dist` (or `mean_spatial_dist_scaled` when `scale = TRUE`).}
#'   \item{value}{Mean value distance from pixels to the supercell center in
#'   value space. Returned as `mean_value_dist` (or `mean_value_dist_scaled`
#'   when `scale = TRUE`).}
#'   \item{combined}{Mean combined distance using `compactness` and `step`
#'   to scale value and spatial distances. Returned as `mean_combined_dist`.}
#'   \item{balance}{Absolute log ratio of scaled value distance to scaled
#'   spatial distance; 0 indicates balance, larger values indicate dominance.}
#' }
#' @seealso [`sc_slic()`], [`sc_metrics_pixels()`], [`sc_metrics_global()`]
#' @export
#' @examples
#' library(supercells)
#' vol = terra::rast(system.file("raster/volcano.tif", package = "supercells"))
#' vol_sc = sc_slic(vol, step = 8, compactness = 7)
#' cl = sc_metrics_supercells(vol, vol_sc)
#' head(cl)
sc_metrics_supercells = function(raster, x, dist_fun = "euclidean", scale = TRUE,
                               metrics = c("spatial", "value", "combined", "balance"),
                               compactness, step) {
  if (any(!metrics %in% c("spatial", "value", "combined", "balance"))) {
    stop("metrics must be one or more of: spatial, value, combined, balance", call. = FALSE)
  }

  prep = .sc_metrics_prep(raster, x, dist_fun, compactness, step)
  x_df = sf::st_drop_geometry(prep$x)
  order_idx = order(x_df[["supercells"]])

  out = sc_metrics_supercells_cpp(prep$clusters, prep$centers_xy, prep$centers_vals, prep$vals,
                                step = prep$step, compactness = prep$compactness,
                                dist_name = prep$dist_name, dist_fun = prep$dist_fun)

  mean_value_dist = out[["mean_value_dist"]]
  mean_spatial_dist = out[["mean_spatial_dist"]]
  mean_combined_dist = out[["mean_combined_dist"]]
  balance = abs(log(out[["balance"]]))

  if (isTRUE(scale)) {
    mean_value_dist = mean_value_dist / prep$compactness
    mean_spatial_dist = mean_spatial_dist / prep$step
  }

  x_ordered = prep$x[order_idx, , drop = FALSE]
  x_keep = x_ordered[, "supercells", drop = FALSE]

  metric_values = list(
    spatial = mean_spatial_dist,
    value = mean_value_dist,
    combined = mean_combined_dist,
    balance = balance
  )
  out_metrics = data.frame(metric_values[metrics])
  name_map = c(
    spatial = if (isTRUE(scale)) "mean_spatial_dist_scaled" else "mean_spatial_dist",
    value = if (isTRUE(scale)) "mean_value_dist_scaled" else "mean_value_dist",
    combined = "mean_combined_dist",
    balance = "balance"
  )
  names(out_metrics) = unname(name_map[metrics])
  return(cbind(x_keep, out_metrics))
}
