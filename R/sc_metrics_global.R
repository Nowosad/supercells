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
#' @param scale Logical. If `TRUE`, scales spatial and value distances; output
#' columns are named with the `_scaled` suffix.
#' @param metrics Character vector of metric ideas to return. Options:
#' `"spatial"`, `"value"`, `"combined"`, `"balance"`. Default:
#' `c("spatial", "value", "combined", "balance")`.
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
#'   \item{mean_value_dist}{Mean per-supercell value distance from cells to their
#'   supercell centers, averaged across supercells. Returned as `mean_value_dist`
#'   (or `mean_value_dist_scaled` when `scale = TRUE`).}
#'   \item{mean_spatial_dist}{Mean per-supercell spatial distance from cells to
#'   their supercell centers, averaged across supercells; units are grid cells
#'   (row/column index distance), not map units. Returned as `mean_spatial_dist`
#'   (or `mean_spatial_dist_scaled` when `scale = TRUE`).}
#'   \item{mean_combined_dist}{Mean per-supercell combined distance, computed from
#'   value and spatial distances using `compactness` and `step`, averaged across
#'   supercells. Returned as `mean_combined_dist`.}
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
                             metrics = c("spatial", "value", "combined", "balance"),
                             compactness, step) {
  if (any(!metrics %in% c("spatial", "value", "combined", "balance"))) {
    stop("metrics must be one or more of: spatial, value, combined, balance", call. = FALSE)
  }

  prep = .sc_metrics_prep(raster, x, dist_fun, compactness, step)
  out = sc_metrics_global_cpp(prep$clusters, prep$centers_xy, prep$centers_vals, prep$vals,
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
  results = cbind(
    data.frame(step = prep$step, compactness = prep$compactness, n_supercells = out[["n_supercells"]]),
    out_metrics
  )
  return(results)
}
