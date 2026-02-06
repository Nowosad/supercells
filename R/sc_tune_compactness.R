#' Estimate compactness from a light SLIC run
#'
#' Runs a short SLIC segmentation (default `iter = 1`) and uses cell-level
#' distances to estimate a compactness value where value and spatial distances
#' are balanced for the chosen `step`.
#' The global estimate uses a pixel-weighted median over the sampled cells,
#' while the local estimate uses a median of per-center mean distances.
#'
#' @param raster A `SpatRaster`.
#' @param step The distance (number of cells) between initial centers (alternative is `k`).
#' @param step_unit Units for `step`. Use "cells" for pixel units or "map" for map units.
#' @param compactness Starting compactness used for the initial short run.
#' @param metrics Which compactness metric to return: `"global"` or `"local"`.
#' Default: `"global"`.
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
#' @return A one-row data frame with columns `step`, `metric`, and `compactness`.
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
sc_tune_compactness = function(raster, step = NULL, step_unit = "cells", compactness = 1,
                                        metrics = "global",
                                        dist_fun = "euclidean", avg_fun = "mean",
                                        clean = TRUE, minarea, iter = 1,
                                        k = NULL, centers = NULL,
                                        sample_size = 10000, value_scale = "auto") {
  pts = sc_slic_points(raster, step = step, compactness = compactness,
                       dist_fun = dist_fun, avg_fun = avg_fun,
                       clean = clean, minarea = minarea, iter = iter,
                       step_unit = step_unit, k = k, centers = centers,
                       outcomes = c("supercells", "coordinates", "values"),
                       chunks = FALSE)

  step_used = attr(pts, "step")
  if (is.null(step_used)) {
    step_used = step
  }

  if (!is.character(metrics) || length(metrics) != 1 || is.na(metrics) ||
      !(metrics %in% c("global", "local"))) {
    stop("metrics must be 'global' or 'local'", call. = FALSE)
  }

  if (identical(value_scale, "auto")) {
    value_scale = sqrt(terra::nlyr(raster))
  } else if (!is.numeric(value_scale) || length(value_scale) != 1 ||
             is.na(value_scale) || value_scale <= 0) {
    stop("value_scale must be a single positive number or 'auto'", call. = FALSE)
  }

  if (identical(metrics, "global")) {
    pix_metrics = sc_metrics_pixels(raster, pts, dist_fun = dist_fun,
                                    compactness = compactness, step = step_used,
                                    scale = FALSE, metrics = c("spatial", "value"))

    med = terra::global(pix_metrics, stats::median, na.rm = TRUE, maxcell = sample_size)
    spatial_dist_median = med[1, 1]
    value_dist_median = med[2, 1]
    value_dist_median = value_dist_median / value_scale
    compactness_value = value_dist_median * step_used / spatial_dist_median

    result = data.frame(step = step_used, metric = "global", compactness = compactness_value)
  }

  if (identical(metrics, "local")) {
    prep = .sc_metrics_prep(raster, pts, dist_fun, compactness, step_used,
                            include = c("centers", "vals", "dist", "raster"))
    mean_value_dist = sc_metrics_local_mean_cpp(
      prep$centers_xy, prep$centers_vals, prep$vals,
      rows = dim(prep$raster)[1], cols = dim(prep$raster)[2],
      step = step_used,
      dist_name = prep$dist_name,
      dist_fun = prep$dist_fun
    )
    mean_value_dist = mean_value_dist / value_scale
    compactness_value = stats::median(mean_value_dist, na.rm = TRUE)

    result = data.frame(step = step_used, metric = "local", compactness = compactness_value)
  }

  return(result)
}
