#' Estimate compactness from a light SLIC run
#'
#' Runs a short SLIC segmentation (default `iter = 1`) and uses cell-level
#' distances to estimate a compactness value where value and spatial distances
#' are balanced for the chosen `step`.
#' Summaries are cell-weighted medians over the (sampled) cells.
#'
#' @param raster A `SpatRaster`.
#' @param step The distance (number of cells) between initial centers (alternative is `k`).
#' @param compactness Starting compactness for the pilot run.
#' @param value_scale Optional scale factor applied to the median value distance
#' before computing compactness. Use `"auto"` to divide by `sqrt(nlyr(raster))`
#' (useful for high-dimensional embeddings). Default: `"auto"`.
#' @param dist_fun A distance function name or a custom function. Supported names:
#' "euclidean", "jsd", "dtw", "dtw2d", or any method from `philentropy::getDistMethods()`.
#' A custom function must accept two numeric vectors and return a single numeric value.
#' @param avg_fun An averaging function name or custom function used to summarize
#' values within each supercell. Supported names: "mean" and "median". A custom
#' function must accept a numeric vector and return a single numeric value.
#' @param clean Should connectivity of the supercells be enforced?
#' @param minarea Minimal size of a supercell (in cells).
#' @param iter Number of SLIC iterations for the pilot run.
#' @param k The number of supercells desired (alternative to `step`).
#' @param centers Optional sf object of custom centers. Requires `step`.
#' @param sample_size Optional limit on the number of pixels used for the summary
#' (passed to `terra::global()` as `maxcell`).
#'
#' @return A one-row data frame with columns `step`, `median_value_dist`,
#' `median_spatial_dist`, and `compactness`.
#'
#' @seealso [`sc_slic()`]
#'
#' @examples
#' library(terra)
#' vol = rast(system.file("raster/volcano.tif", package = "supercells"))
#' tune = sc_tune_compactness(vol, step = 8)
#' tune$compactness
#'
#' @export
sc_tune_compactness = function(raster, step = NULL, compactness = 1,
                                        dist_fun = "euclidean", avg_fun = "mean",
                                        clean = TRUE, minarea, iter = 1,
                                        k = NULL, centers = NULL,
                                        sample_size = 10000, value_scale = "auto") {
  pts = sc_slic_points(raster, step = step, compactness = compactness,
                       dist_fun = dist_fun, avg_fun = avg_fun,
                       clean = clean, minarea = minarea, iter = iter,
                       k = k, centers = centers, metadata = TRUE,
                       chunks = FALSE)

  step_used = if (is.null(step)) attr(pts, "step") else step

  metrics = sc_metrics_pixels(raster, pts, dist_fun = dist_fun,
                              compactness = compactness, step = step_used,
                              scale = FALSE, metrics = c("spatial", "value"))

  if (missing(sample_size)) {
    sample_size = Inf
  }
  med = terra::global(metrics, stats::median, na.rm = TRUE, maxcell = sample_size)
  median_spatial_dist = med[1, 1]
  median_value_dist = med[2, 1]

  if (is.character(value_scale)) {
    if (identical(value_scale, "auto")) {
      value_scale = sqrt(terra::nlyr(raster))
    } else {
      stop("value_scale must be numeric or 'auto'", call. = FALSE)
    }
  }
  if (!is.numeric(value_scale) || length(value_scale) != 1 || is.na(value_scale) || value_scale <= 0) {
    stop("value_scale must be a single positive number or 'auto'", call. = FALSE)
  }
  median_value_dist = median_value_dist / value_scale

  compactness = median_value_dist * step_used / median_spatial_dist

  data.frame(
    step = step_used,
    median_value_dist = median_value_dist,
    median_spatial_dist = median_spatial_dist,
    compactness = compactness
  )
}
