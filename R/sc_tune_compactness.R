#' Estimate a compactness value
#'
#' Estimates a compactness value for a chosen raster scale.
#' The current implementation supports one tuning metric, `"local_variability"`,
#' which estimates compactness directly from local value variability.
#'
#' @param raster A `SpatRaster`.
#' @param step Initial center spacing (alternative is `k`).
#' Provide a plain numeric value for cell units, or use [use_meters()] for
#' map-distance steps in meters (automatically converted to cells using raster resolution).
#' @param dist_fun A distance function name or a custom function. Supported names:
#' "euclidean", "jsd", "dtw", "dtw2d", or any method from `philentropy::getDistMethods()`.
#' A custom function must accept two numeric vectors and return a single numeric value.
#' @param metric Which compactness metric to return. Currently only
#' `"local_variability"` is supported. The argument is kept for future additions.
#' For `"local_variability"`, `compactness = median(local_mean(d_value)) / dim_scale`,
#' where local means are computed in windows around initial centers and
#' `dim_scale` is inferred from the number of raster layers and `dist_fun`.
#' This keeps compactness adjustable to dimensionality without requiring a
#' separate user-facing scaling argument.
#' @param k The number of supercells desired (alternative to `step`).
#' @param centers Optional sf object of custom initial centers. Requires `step`.
#'
#' @return A one-row data frame with columns `step` and `compactness`.
#'
#' @seealso [`sc_slic()`], [use_meters()], [use_adaptive()]
#'
#' @examples
#' library(terra)
#' vol = rast(system.file("raster/volcano.tif", package = "supercells"))
#' tune = sc_tune_compactness(vol, step = 8)
#' tune$compactness
#'
#' @export
sc_tune_compactness = function(raster, step = NULL, dist_fun = "euclidean",
                                        metric = "local_variability", k = NULL,
                                        centers = NULL) {
  if (!is.character(metric) || length(metric) != 1 || is.na(metric) ||
      !identical(metric, "local_variability")) {
    stop("metric must be 'local_variability'", call. = FALSE)
  }
  dist_fun_out = if (is.character(dist_fun)) as.character(dist_fun[[1]]) else "custom"

  dim_scale = .sc_tune_local_variability_scale(terra::nlyr(raster), dist_fun_out)

  tune_stats = .sc_tune_grid_window_variability(raster, step, dist_fun, k, centers)
  step_used = unlist(tune_stats$step, use.names = FALSE)
  compactness_value = tune_stats$compactness / dim_scale

  return(data.frame(step = step_used, metric = "local_variability", dist_fun = dist_fun_out,
                    compactness = compactness_value))
}

.sc_tune_local_variability_scale = function(n_bands, dist_fun_name) {
  if (!is.character(dist_fun_name) || length(dist_fun_name) != 1 || is.na(dist_fun_name)) {
    return(1)
  }

  if (identical(dist_fun_name, "euclidean")) {
    return(sqrt(n_bands))
  }

  if (dist_fun_name %in% c("manhattan", "dtw", "dtw2d")) {
    return(n_bands)
  }

  if (dist_fun_name %in% c("jsd")) {
    return(1)
  }

  1
}

.sc_tune_grid_window_variability = function(raster, step, dist_fun, k, centers) {
  pts = sc_slic_points(
    raster,
    step = step,
    compactness = 1,
    dist_fun = dist_fun,
    iter = 0,
    k = k,
    centers = centers,
    outcomes = c("supercells", "coordinates", "values"),
    chunks = FALSE
  )

  step_used = attr(pts, "step")
  prep = .sc_metrics_prep(raster, pts, dist_fun, compactness = 1, step = step_used,
                          include = c("centers", "vals", "dist", "raster"))

  mean_value_dist = sc_metrics_local_mean_cpp(
    prep$centers_xy, prep$centers_vals, prep$vals,
    rows = dim(prep$raster)[1], cols = dim(prep$raster)[2],
    step = prep$step,
    dist_name = prep$dist_name,
    dist_fun = prep$dist_fun
  )

  list(
    step = attr(pts, "step"),
    compactness = stats::median(mean_value_dist, na.rm = TRUE)
  )
}
