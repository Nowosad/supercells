# .sc_metrics_prep: normalize inputs and parameters for metrics functions
# Inputs: raster and supercells; outputs include prepared matrices and metadata
# Handles missing metadata by deriving centers and ids from geometry
.sc_metrics_prep = function(raster, x, dist_fun, compactness, step,
                            include = c("clusters", "centers", "vals", "dist", "raster")) {

  valid = c("clusters", "centers", "vals", "dist", "raster")
  if (!all(include %in% valid)) {
    stop("include must be one or more of: ", paste(valid, collapse = ", "), call. = FALSE)
  }

  # prepare arguments
  raster = .sc_util_prep_raster(raster)

  if (!inherits(x, "sf")) {
    stop("The 'x' argument must be an sf object returned by sc_slic()", call. = FALSE)
  }
  adaptive_compactness = FALSE
  if (missing(compactness)) {
    method = attr(x, "method")
    adaptive_compactness = isTRUE(identical(method, "slic0"))
    compactness = attr(x, "compactness")
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
  result = list(
    x = x_work,
    step = step,
    compactness = compactness,
    adaptive_compactness = adaptive_compactness
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
    cluster_rast = terra::rasterize(terra::vect(x_work), raster, field = "supercells")
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
