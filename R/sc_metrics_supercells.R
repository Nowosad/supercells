#' Supercell-level metrics
#'
#' Computes per-supercell distance diagnostics
#'
#' @inheritParams sc_metrics_pixels
#' @param metrics Character vector of metrics to return. Options:
#' `"mean_spatial_dist"`, `"mean_value_dist"`, `"mean_combined_dist"`,
#' `"balance"`. Default:
#' `c("mean_spatial_dist", "mean_value_dist", "mean_combined_dist", "balance")`.
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
#'   \item{mean_value_dist}{Mean value distance from pixels to the supercell
#'   center in value space.}
#'   \item{mean_spatial_dist}{Mean spatial distance from pixels to the
#'   supercell center in grid-cell units (row/column index distance).}
#'   \item{mean_combined_dist}{Mean combined distance using `compactness`
#'   and `step` to scale value and spatial distances.}
#'   \item{balance}{Absolute log ratio of scaled value distance to
#'   scaled spatial distance; 0 indicates balance, larger values indicate dominance.}
#' }
#' When `scale = TRUE`, `mean_spatial_dist` and `mean_value_dist` are returned as
#' `mean_spatial_dist_scaled` and `mean_value_dist_scaled`.
#' @seealso [`sc_slic()`], [`sc_metrics_pixels()`], [`sc_metrics_global()`]
#' @export
#' @examples
#' library(supercells)
#' vol = terra::rast(system.file("raster/volcano.tif", package = "supercells"))
#' vol_sc = sc_slic(vol, step = 8, compactness = 7)
#' cl = sc_metrics_supercells(vol, vol_sc)
#' head(cl)
sc_metrics_supercells = function(raster, x, dist_fun = "euclidean", scale = TRUE,
                               metrics = c("mean_spatial_dist", "mean_value_dist",
                                           "mean_combined_dist", "balance"),
                               compactness, step) {
  prep = .sc_metrics_prep(raster, x, dist_fun, compactness, step)
  x_df = sf::st_drop_geometry(prep$x)
  order_idx = order(x_df[["supercells"]])

  out = sc_metrics_supercells_cpp(prep$clusters, prep$centers_xy, prep$centers_vals, prep$vals,
                                step = prep$step, compactness = prep$compactness,
                                dist_name = prep$dist_name, dist_fun = prep$dist_fun)

  metrics_df = data.frame(
    mean_value_dist = out[["mean_value_dist"]],
    mean_spatial_dist = out[["mean_spatial_dist"]],
    mean_combined_dist = out[["mean_combined_dist"]],
    balance = abs(log(out[["balance"]]))
  )
  if (isTRUE(scale)) {
    metrics_df$mean_value_dist = metrics_df$mean_value_dist / prep$compactness
    metrics_df$mean_spatial_dist = metrics_df$mean_spatial_dist / prep$step
  }
  if (any(!metrics %in% names(metrics_df))) {
    bad = metrics[!metrics %in% names(metrics_df)]
    stop(sprintf("Unknown metrics: %s", paste(bad, collapse = ", ")), call. = FALSE)
  }

  x_ordered = prep$x[order_idx, , drop = FALSE]
  x_keep = x_ordered[, "supercells", drop = FALSE]
  out_metrics = metrics_df[, metrics, drop = FALSE]
  if (isTRUE(scale)) {
    names(out_metrics) = sub("^mean_spatial_dist$", "mean_spatial_dist_scaled", names(out_metrics))
    names(out_metrics) = sub("^mean_value_dist$", "mean_value_dist_scaled", names(out_metrics))
  }
  return(cbind(x_keep, out_metrics))
}
