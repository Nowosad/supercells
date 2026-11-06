#' Estimate a compactness value
#'
#' Estimates a compactness value for a chosen raster scale.
#' The current implementation supports one tuning metric, `"local_variability"`,
#' which estimates compactness directly from local value variability.
#'
#' @param x A `SpatRaster` (terra) or `stars` object.
#' @param step Initial center spacing (alternative is `k`).
#' Provide a plain numeric value for cell units, or use [use_meters()] for
#' map-distance steps in meters (automatically converted to cells using raster resolution).
#' @param dist_fun A supported distance name from [sc_list_distances()] or a
#' custom function that takes two numeric vectors and returns one finite,
#' nonnegative value. Smaller values must mean greater similarity.
#' @param metric Which compactness metric to return. Currently only
#' `"local_variability"` is supported. The argument is kept for future additions.
#' For `"local_variability"`, `compactness = median(local_mean(d_value)) / dim_scale`,
#' where local means are computed in windows around initial centers and
#' `dim_scale` is inferred from the number of raster layers and `dist_fun`.
#' Root-sum-square measures use `sqrt(n_layers)`, additive and path-sum
#' measures use `n_layers`, and other named or custom measures use `1`.
#' This keeps compactness adjustable to dimensionality without requiring a
#' separate user-facing scaling argument.
#' @param stat The summary statistic to use for the local variability metric. Default is `median`. If it is `NULL`, all local variability values will be returned without summarization.
#' @param k The number of supercells desired (alternative to `step`).
#' @param centers Optional sf object of custom initial centers. Requires `step`.
#' @param ... Additional arguments for the compactness metric. For example, `dim_scale` for `"local_variability"` if you want to override the default scaling.
#'
#' @return A data frame with columns `step`, `metric`, `dist_fun`, and
#' `compactness`. With the default `stat = median`, it contains one row. With
#' `stat = NULL`, it contains one row per initial center, with `step`, `metric`,
#' and `dist_fun` repeated for each row.
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
sc_tune_compactness = function(x, step = NULL, dist_fun = "euclidean",
                                metric = "local_variability", stat = stats::median,
                                k = NULL, centers = NULL, ...) {
  if (!is.character(metric) || length(metric) != 1 || is.na(metric) ||
      !identical(metric, "local_variability")) {
    stop("metric must be 'local_variability'", call. = FALSE)
  }
  x = .sc_util_prep_raster(x)
  dist_fun_out = if (is.character(dist_fun)) as.character(dist_fun[[1]]) else "custom"

  tune_stats = .sc_tune_grid_window_variability(x, step, dist_fun, k, centers)
  step_used = tune_stats$step

  if (!is.null(list(...)$dim_scale)) {
    dim_scale = list(...)$dim_scale
  } else {
    dim_scale = .sc_tune_local_variability_scale(terra::nlyr(x), dist_fun_out)
  }
  compactness_value = tune_stats$compactness / dim_scale

  if (!is.null(stat)) {
    compactness_value = compactness_value[!is.na(compactness_value)]
    compactness_value = stat(compactness_value)
  }

  n_out = length(compactness_value)
  return(data.frame(
    step = rep(step_used, length.out = n_out),
    metric = rep("local_variability", length.out = n_out),
    dist_fun = rep(dist_fun_out, length.out = n_out),
    compactness = compactness_value
  ))
}

.sc_tune_local_variability_scale = function(n_bands, dist_fun_name) {
  if (!is.character(dist_fun_name) || length(dist_fun_name) != 1 || is.na(dist_fun_name)) {
    return(1)
  }

  tuning_scale = .sc_distance_registry$tuning_scale[
    match(dist_fun_name, .sc_distance_registry$name)
  ]
  if (is.na(tuning_scale)) return(1)

  switch(
    tuning_scale,
    none = 1,
    sqrt_layers = sqrt(n_bands),
    layers = n_bands,
    stop("Unknown compactness tuning scale", call. = FALSE)
  )
}

.sc_tune_grid_window_variability = function(x, step, dist_fun, k, centers) {
  pts = sc_slic_points(
    x,
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
  prep = .sc_metrics_prep(x, pts, dist_fun, compactness = 1, step = step_used,
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
    compactness = mean_value_dist
  )
}
