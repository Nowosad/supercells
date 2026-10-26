# .sc_metrics_prep: normalize inputs and parameters for metrics functions
# Inputs: raster and supercells; outputs include prepared matrices and metadata
# Handles missing metadata by deriving centers and ids from geometry

.sc_metrics_resolve_dist_fun = function(sc, dist_fun) {
  if (!is.null(dist_fun)) {
    return(dist_fun)
  }
  dist_fun = attr(sc, "dist_fun")
  if (is.null(dist_fun)) {
    stop("The 'dist_fun' argument is required when it is not stored in 'sc'", call. = FALSE)
  }
  dist_fun
}

.sc_metrics_scale_summary = function(value_dist, spatial_dist, out, prep, scale) {
  if (!isTRUE(scale)) {
    return(list(value_dist = value_dist, spatial_dist = spatial_dist))
  }
  if (isTRUE(prep$adaptive_compactness)) {
    value_dist = out[["mean_value_dist_scaled"]]
  } else {
    value_dist = value_dist / prep$compactness
  }
  spatial_dist = spatial_dist / prep$step_scale
  list(value_dist = value_dist, spatial_dist = spatial_dist)
}

.sc_metrics_prep = function(x, sc, dist_fun, compactness, step,
                            include = c("clusters", "centers", "vals", "dist", "raster")) {

  valid = c("clusters", "centers", "vals", "dist", "raster")
  if (!all(include %in% valid)) {
    stop("include must be one or more of: ", paste(valid, collapse = ", "), call. = FALSE)
  }

  # prepare arguments
  raster = .sc_util_prep_raster(x)

  if (!inherits(sc, "sf")) {
    stop("The 'sc' argument must be an sf object returned by sc_slic()", call. = FALSE)
  }

  if (missing(compactness)) {
    compactness = attr(sc, "compactness")
  }
  compactness_method = attr(sc, "compactness_method")
  if (missing(step)) {
    step = attr(sc, "step")
  }
  if (is.null(step)) {
    stop("Both 'compactness' and 'step' are required", call. = FALSE)
  }

  if (!is.null(compactness_method)) {
    if (!is.character(compactness_method) || length(compactness_method) != 1 || is.na(compactness_method) ||
        !(compactness_method %in% c("constant", "local_max"))) {
      stop("The 'compactness_method' attribute must be 'constant' or 'local_max'", call. = FALSE)
    }
  }

  if (identical(compactness_method, "local_max")) {
    compactness_prep = list(value = 0, adaptive = TRUE, compactness_method = "local_max")
  } else {
    compactness_prep = .sc_util_prep_compactness(compactness)
    if (!is.null(compactness_method) && compactness_prep$compactness_method != compactness_method) {
      stop("The provided compactness method conflicts with 'compactness_method' attribute", call. = FALSE)
    }
  }
  step_prep = .sc_util_step_to_cells(raster, step)
  step = step_prep$step

  # prepare data, including handling missing metadata
  sc_work = sc
  x_df = sf::st_drop_geometry(sc_work)
  if (!("supercells" %in% names(x_df))) {
    sc_work[["supercells"]] = seq_len(nrow(sc_work))
    x_df = sf::st_drop_geometry(sc_work)
  }
  if (!all(c("x", "y") %in% names(x_df))) {
    old_s2 = sf::sf_use_s2()
    on.exit(suppressMessages(sf::sf_use_s2(old_s2)), add = TRUE)
    suppressMessages(sf::sf_use_s2(FALSE))
    centers = sf::st_centroid(sf::st_geometry(sc_work))
    coords = sf::st_coordinates(centers)
    x_df[["x"]] = coords[, 1]
    x_df[["y"]] = coords[, 2]
  }
  val_cols = names(raster)
  if (!all(val_cols %in% names(x_df))) {
    stop("The value columns in 'sc' must match the raster layer names", call. = FALSE)
  }
  x_df = x_df[order(x_df[["supercells"]]), , drop = FALSE]

  # prepare matrices for C++ function
  spatial_scale = step_prep$spatial_scale
  step_scale = step_prep$step_scale

  result = list(
    sc = sc_work,
    step = step,
    step_meta = step_prep$step_meta,
    compactness = compactness_prep$value,
    adaptive_compactness = compactness_prep$adaptive,
    compactness_method = compactness_prep$compactness_method,
    spatial_scale = spatial_scale,
    step_scale = step_scale
  )

  if ("centers" %in% include) {
    raster_ext = as.vector(terra::ext(raster))
    raster_res = terra::res(raster)
    center_x = (x_df[["x"]] - raster_ext[[1]]) / raster_res[[1]] - 0.5
    center_y = (raster_ext[[4]] - x_df[["y"]]) / raster_res[[2]] - 0.5
    tolerance = sqrt(.Machine$double.eps)
    outside = !is.finite(center_x) | !is.finite(center_y) |
      center_x < -tolerance | center_x > terra::ncol(raster) - 1 + tolerance |
      center_y < -tolerance | center_y > terra::nrow(raster) - 1 + tolerance
    if (any(outside)) {
      stop("Some centers fall outside the raster extent", call. = FALSE)
    }
    center_x = pmin(pmax(center_x, 0), terra::ncol(raster) - 1)
    center_y = pmin(pmax(center_y, 0), terra::nrow(raster) - 1)
    centers_xy = cbind(center_x, center_y)
    storage.mode(centers_xy) = "double"
    centers_vals = as.matrix(x_df[, val_cols, drop = FALSE])
    storage.mode(centers_vals) = "double"
    result$centers_xy = centers_xy
    result$centers_vals = centers_vals
  }

  if ("clusters" %in% include) {
    internal_ids = match(sc_work[["supercells"]], x_df[["supercells"]])
    cluster_rast = terra::rasterize(terra::vect(sc_work), raster, field = internal_ids)
    clusters = terra::as.matrix(cluster_rast, wide = TRUE)
    clusters = ifelse(is.na(clusters), -1L, clusters - 1L)
    storage.mode(clusters) = "integer"
    result$clusters = clusters
  }

  if ("vals" %in% include) {
    vals = terra::values(raster, mat = TRUE, na.rm = FALSE)
    storage.mode(vals) = "double"
    result$vals = vals
  }

  if ("dist" %in% include) {
    dist_prep = .sc_util_prep_dist_fun(dist_fun)
    result$dist_name = dist_prep$dist_name
    result$dist_fun = dist_prep$dist_fun
  }

  if ("raster" %in% include) {
    result$raster = raster
  }

  return(result)
}
