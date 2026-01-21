#' Global supercells metrics
#'
#' Computes global distance diagnostics for supercells.
#'
#' @param x An sf object returned by [supercells()].
#' @param raster The input SpatRaster used to create `x`.
#' @param dist_fun A distance function name or function, as in [supercells()].
#' @param compactness A compactness value used for the supercells.
#' @param step A step value used for the supercells.
#'
#' @return A data.frame with a single row of global metrics.
#' @export
sc_metrics_global = function(x, raster, dist_fun = "euclidean", compactness, step) {
  if (!inherits(raster, "SpatRaster")) {
    raster = terra::rast(raster)
  }
  if (!inherits(x, "sf")) {
    stop("The 'x' argument must be an sf object returned by supercells()", call. = FALSE)
  }
  if (missing(compactness) || missing(step)) {
    stop("Both 'compactness' and 'step' are required", call. = FALSE)
  }

  x_df = sf::st_drop_geometry(x)
  if (!all(c("supercells", "x", "y") %in% names(x_df))) {
    stop("The 'x' object must contain 'supercells', 'x', and 'y' columns", call. = FALSE)
  }
  val_cols = setdiff(names(x_df), c("supercells", "x", "y"))
  if (length(val_cols) == 0) {
    stop("No value columns found in 'x'", call. = FALSE)
  }
  x_df = x_df[order(x_df[["supercells"]]), , drop = FALSE]

  center_x = terra::colFromX(raster, x_df[["x"]])
  center_y = terra::rowFromY(raster, x_df[["y"]])
  if (any(is.na(center_x)) || any(is.na(center_y))) {
    stop("Some centers fall outside the raster extent", call. = FALSE)
  }
  centers_xy = cbind(center_x, center_y)
  storage.mode(centers_xy) = "double"
  centers_vals = as.matrix(x_df[, val_cols, drop = FALSE])
  storage.mode(centers_vals) = "double"

  cluster_rast = terra::rasterize(terra::vect(x), raster, field = "supercells")
  clusters = terra::as.matrix(cluster_rast, wide = TRUE)
  clusters = ifelse(is.na(clusters), -1L, clusters - 1L)
  storage.mode(clusters) = "integer"

  vals = as.matrix(terra::as.data.frame(raster, cells = FALSE, na.rm = FALSE))
  storage.mode(vals) = "double"

  if (is.character(dist_fun)) {
    if (!(dist_fun %in% c("euclidean", "jsd", "dtw", "dtw2d", philentropy::getDistMethods()))) {
      stop("The provided distance function ('dist_fun') does not exist!", call. = FALSE)
    }
    dist_name = dist_fun
    dist_fun = function() ""
  } else {
    dist_name = ""
  }

  out = sc_metrics_global_cpp(clusters, centers_xy, centers_vals, vals,
                              step = step, compactness = compactness,
                              dist_name = dist_name, dist_fun = dist_fun)

  data.frame(
    step = step,
    compactness = compactness,
    n_supercells = out[["n_supercells"]],
    mean_value_dist = out[["mean_value_dist"]],
    mean_spatial_dist = out[["mean_spatial_dist"]],
    mean_combined_dist = out[["mean_combined_dist"]],
    compactness_ratio_mean = out[["compactness_ratio_mean"]],
    mean_value_dist_w = out[["mean_value_dist_w"]],
    mean_spatial_dist_w = out[["mean_spatial_dist_w"]],
    mean_combined_dist_w = out[["mean_combined_dist_w"]]
  )
}
