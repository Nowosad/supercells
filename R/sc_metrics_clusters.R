#' Cluster-level supercells metrics
#'
#' Computes per-cluster distance diagnostics
#'
#' @inheritParams sc_metrics_pixels
#' @details
#' If `x` lacks `supercells`, `x`, or `y` columns, they are derived from geometry
#' and row order, which may differ from the original centers
#' @return An sf object with one row per supercell
#' @export
#' @examples
#' library(supercells)
#' vol = terra::rast(system.file("raster/volcano.tif", package = "supercells"))
#' vol_sc = sc_slic(vol, step = 8, compactness = 1, metadata = TRUE)
#' cl = sc_metrics_clusters(vol, vol_sc)
#' head(cl)
sc_metrics_clusters = function(raster, x, dist_fun = "euclidean", compactness, step) {
  prep = .sc_metrics_prep(raster, x, dist_fun, compactness, step)
  x_df = sf::st_drop_geometry(prep$x)
  order_idx = order(x_df[["supercells"]])

  out = sc_metrics_clusters_cpp(prep$clusters, prep$centers_xy, prep$centers_vals, prep$vals,
                                step = prep$step, compactness = prep$compactness,
                                dist_name = prep$dist_name, dist_fun = prep$dist_fun)

  metrics = data.frame(
    mean_value_dist = out[["mean_value_dist"]],
    mean_spatial_dist = out[["mean_spatial_dist"]],
    mean_combined_dist = out[["mean_combined_dist"]],
    compactness_ratio = out[["compactness_ratio"]]
  )

  x_ordered = prep$x[order_idx, , drop = FALSE]
  x_keep = x_ordered[, "supercells", drop = FALSE]
  return(cbind(x_keep, metrics))
}
