# .sc_metrics_prep: normalize inputs and parameters for metrics functions
# Inputs: raster and supercells; outputs include prepared matrices and metadata
# Handles missing metadata by deriving centers and ids from geometry
.sc_metrics_prep = function(raster, x, dist_fun, compactness, step) {

  # prepare arguments
  raster = .sc_util_prep_raster(raster)

  if (!inherits(x, "sf")) {
    stop("The 'x' argument must be an sf object returned by sc_slic()", call. = FALSE)
  }
  adaptive_compactness = FALSE
  if (missing(compactness)) {
    compactness = attr(x, "compactness")
    method = attr(x, "method")
    adaptive_compactness = isTRUE(identical(method, "slic0"))
  } else if (is.character(compactness) && length(compactness) == 1 && !is.na(compactness) && compactness == "auto") {
    adaptive_compactness = TRUE
    compactness = 0
  }
  if (missing(step)) {
    step = attr(x, "step")
  }
  if (is.null(compactness) || is.null(step)) {
    stop("Both 'compactness' and 'step' are required", call. = FALSE)
  }

  # prepare data, including handling missing metadata
  x_work = x
  x_df = sf::st_drop_geometry(x_work)
  if (!("supercells" %in% names(x_df))) {
    x_work[["supercells"]] = seq_len(nrow(x_work))
    x_df = sf::st_drop_geometry(x_work)
  }
  if (!all(c("x", "y") %in% names(x_df))) {
    centers = sf::st_centroid(sf::st_geometry(x_work))
    coords = sf::st_coordinates(centers)
    x_df[["x"]] = coords[, 1]
    x_df[["y"]] = coords[, 2]
  }
  val_cols = setdiff(names(x_df), c("supercells", "x", "y"))
  if (length(val_cols) == 0) {
    stop("No value columns found in 'x'", call. = FALSE)
  }
  x_df = x_df[order(x_df[["supercells"]]), , drop = FALSE]

  # prepare matrices for C++ function
  center_x = terra::colFromX(raster, x_df[["x"]])
  center_y = terra::rowFromY(raster, x_df[["y"]])
  if (any(is.na(center_x)) || any(is.na(center_y))) {
    stop("Some centers fall outside the raster extent", call. = FALSE)
  }
  centers_xy = cbind(center_x, center_y)
  storage.mode(centers_xy) = "double"

  centers_vals = as.matrix(x_df[, val_cols, drop = FALSE])
  storage.mode(centers_vals) = "double"

  cluster_rast = terra::rasterize(terra::vect(x_work), raster, field = "supercells")
  clusters = terra::as.matrix(cluster_rast, wide = TRUE)
  clusters = ifelse(is.na(clusters), -1L, clusters - 1L)
  storage.mode(clusters) = "integer"

  vals = as.matrix(terra::as.data.frame(raster, cells = FALSE, na.rm = FALSE))
  storage.mode(vals) = "double"

  dist_prep = .sc_util_prep_dist_fun(dist_fun)

  result = list(
    raster = raster,
    x = x_work,
    clusters = clusters,
    centers_xy = centers_xy,
    centers_vals = centers_vals,
    vals = vals,
    step = step,
    compactness = compactness,
    adaptive_compactness = adaptive_compactness,
    dist_name = dist_prep$dist_name,
    dist_fun = dist_prep$dist_fun
  )

  return(result)
}
