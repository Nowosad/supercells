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
  if (missing(step)) {
    step = attr(sc, "step")
  }
  if (is.null(compactness) || is.null(step)) {
    stop("Both 'compactness' and 'step' are required", call. = FALSE)
  }

  compactness_input = compactness
  adaptive_compactness = FALSE
  if (is.character(compactness)) {
    if (length(compactness) != 1 || is.na(compactness) || compactness != "auto") {
      stop("The 'compactness' argument must be numeric or 'auto'", call. = FALSE)
    }
    adaptive_compactness = TRUE
    compactness = 0
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
  val_cols = setdiff(names(x_df), c("supercells", "x", "y"))
  if (length(val_cols) == 0) {
    stop("No value columns found in 'sc'", call. = FALSE)
  }
  x_df = x_df[order(x_df[["supercells"]]), , drop = FALSE]

  # prepare matrices for C++ function
  spatial_scale = step_prep$spatial_scale
  step_scale = step_prep$step_scale

  result = list(
    sc = sc_work,
    step = step,
    step_meta = step_prep$step_meta,
    compactness = compactness,
    compactness_input = compactness_input,
    adaptive_compactness = adaptive_compactness,
    spatial_scale = spatial_scale,
    step_scale = step_scale
  )

  if ("centers" %in% include) {
    center_x = terra::colFromX(raster, x_df[["x"]])
    center_y = terra::rowFromY(raster, x_df[["y"]])
    if (any(is.na(center_x)) || any(is.na(center_y))) {
      stop("Some centers fall outside the raster extent", call. = FALSE)
    }
    centers_xy = cbind(center_x, center_y)
    storage.mode(centers_xy) = "double"
    centers_vals = as.matrix(x_df[, val_cols, drop = FALSE])
    storage.mode(centers_vals) = "double"
    result$centers_xy = centers_xy
    result$centers_vals = centers_vals
  }

  if ("clusters" %in% include) {
    cluster_rast = terra::rasterize(terra::vect(sc_work), raster, field = "supercells")
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
