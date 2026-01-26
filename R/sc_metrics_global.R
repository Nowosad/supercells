#' Global supercells metrics
#'
#' Computes global distance diagnostics for supercells
#'
#' @details
#' Requires `x` with metadata columns (`supercells`, `x`, `y`)
#' If they are missing, they are derived from geometry and row order
#' Set `metadata = TRUE` when calling `sc_slic()` or `supercells()`
#' Metrics are averaged across supercells (each supercell has equal weight).
#'
#' @inheritParams sc_metrics_pixels
#' @param metrics Character vector of metrics to return. Options:
#' `"mean_spatial_dist"`, `"mean_value_dist"`, `"mean_combined_dist"`, `"balance"`.
#' Default: `c("mean_spatial_dist", "mean_value_dist", "mean_combined_dist", "balance")`.
#' @return A data.frame with a single row of global metrics and columns:
#' Interpretation:
#' \describe{
#'   \item{mean_value_dist}{Lower values indicate more homogeneous supercells.}
#'   \item{mean_spatial_dist}{Lower values indicate more compact supercells.}
#'   \item{mean_combined_dist}{Overall distance; mainly useful for ranking.}
#'   \item{balance}{0 indicates balance between value and spatial terms.}
#' }
#' \describe{
#'   \item{step}{Step size used to generate supercells.}
#'   \item{compactness}{Compactness value used to generate supercells.}
#'   \item{n_supercells}{Number of supercells with at least one non-missing pixel.}
#'   \item{mean_value_dist}{Mean per-supercell value distance from pixels to their
#'   supercell centers, averaged across supercells.}
#'   \item{mean_spatial_dist}{Mean per-supercell spatial distance from pixels to
#'   their supercell centers, averaged across supercells; units are grid cells
#'   (row/column index distance), not map units.}
#'   \item{mean_combined_dist}{Mean per-supercell combined distance, computed from
#'   value and spatial distances using `compactness` and `step`, averaged across
#'   supercells.}
#'   \item{balance}{Mean absolute log ratio of scaled value distance to scaled
#'   spatial distance; 0 indicates balance.}
#' }
#' When `scale = TRUE`, `mean_spatial_dist` and `mean_value_dist` are returned as
#' `mean_spatial_dist_scaled` and `mean_value_dist_scaled`.
#' @seealso [`sc_slic()`], [`sc_metrics_pixels()`], [`sc_metrics_supercells()`]
#' @export
#' @examples
#' library(supercells)
#' vol = terra::rast(system.file("raster/volcano.tif", package = "supercells"))
#' vol_sc = sc_slic(vol, step = 8, compactness = 7)
#' sc_metrics_global(vol, vol_sc)
sc_metrics_global = function(raster, x, dist_fun = "euclidean", scale = TRUE,
                             metrics = c("mean_spatial_dist", "mean_value_dist",
                                         "mean_combined_dist", "balance"),
                             compactness, step) {
  prep = .sc_metrics_prep(raster, x, dist_fun, compactness, step)
  out = sc_metrics_global_cpp(prep$clusters, prep$centers_xy, prep$centers_vals, prep$vals,
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
  out_metrics = metrics_df[, metrics, drop = FALSE]
  if (isTRUE(scale)) {
    names(out_metrics) = sub("^mean_spatial_dist$", "mean_spatial_dist_scaled", names(out_metrics))
    names(out_metrics) = sub("^mean_value_dist$", "mean_value_dist_scaled", names(out_metrics))
  }
  results = cbind(
    data.frame(step = prep$step, compactness = prep$compactness, n_supercells = out[["n_supercells"]]),
    out_metrics
  )
  return(results)
}
