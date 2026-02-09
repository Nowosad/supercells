# normalize input to SpatRaster
.sc_util_prep_raster = function(x) {
  if (inherits(x, "SpatRaster")) {
    return(x)
  }
  if (inherits(x, "stars")) {
    return(terra::rast(x))
  }
  stop("The SpatRaster class is expected as an input", call. = FALSE)
}

# validate and prepare distance function
.sc_util_prep_dist_fun = function(dist_fun) {
  if (is.character(dist_fun)) {
    if (!(dist_fun %in% c("euclidean", "jsd", "dtw", "dtw2d", philentropy::getDistMethods()))) {
      stop("The provided distance function ('dist_fun') does not exist!", call. = FALSE)
    }
    dist_name = dist_fun
    dist_fun = function() ""
  } else {
    dist_name = ""
  }
  return(list(dist_name = dist_name, dist_fun = dist_fun))
}

# validate and prepare averaging function
.sc_util_prep_avg_fun = function(avg_fun) {
  if (is.character(avg_fun)) {
    avg_fun_name = avg_fun
    avg_fun_fun = function() ""
  } else {
    avg_fun_name = ""
    avg_fun_fun = avg_fun
  }
  return(list(avg_fun_name = avg_fun_name, avg_fun_fun = avg_fun_fun))
}

# convert sf centers to raster row and col indices
.sc_util_centers_to_dims = function(x, y){
  y_coords = sf::st_coordinates(sf::st_geometry(y))
  y_col = terra::colFromX(x, y_coords[, 1])
  y_row = terra::rowFromY(x, y_coords[, 2])
  center_dims = cbind(y_col, y_row)
  storage.mode(center_dims) = "integer"
  unique(center_dims)
}

# create seq.int including the final to value
.sc_util_seq_last = function(from, to, by){
  vec = do.call(what = seq.int, args = list(from, to, by))
  if (utils::tail(vec, 1) != to) {
    return(c(vec, to))
  } else {
    return(vec)
  }
}

# convert user step to internal cell units and return scaling metadata
.sc_util_step_to_cells = function(x, step) {
  if (!inherits(step, "units")) {
    return(list(step = step, step_meta = step, spatial_scale = 1, step_scale = step))
  }

  if (!identical(as.character(units::deparse_unit(step)), "m")) {
    stop("A units-based 'step' must use meters ('m')", call. = FALSE)
  }
  if (terra::is.lonlat(x)) {
    stop("A units-based 'step' requires a projected CRS; project the input raster first", call. = FALSE)
  }

  res = terra::res(x)
  if (!isTRUE(all.equal(res[[1]], res[[2]]))) {
    warning("Map-unit step requires square cells; res(x) has different x/y resolution", call. = FALSE)
  }

  crs_units = sf::st_crs(terra::crs(x))$units_gdal
  if (is.null(crs_units) || is.na(crs_units) || !nzchar(crs_units)) {
    stop("The raster CRS has unknown linear units; cannot use units-based 'step'", call. = FALSE)
  }
  crs_units_lc = tolower(trimws(crs_units))
  if (!(crs_units_lc %in% c("meter", "metre", "m"))) {
    stop("A units-based 'step' requires a projected CRS with meter units", call. = FALSE)
  }

  step_m = as.numeric(units::drop_units(step))
  step_cells = max(1, round(step_m / res[[1]]))
  step_meta = units::set_units(step_cells * res[[1]], "m", mode = "standard")
  spatial_scale = res[[1]]

  list(step = step_cells, step_meta = step_meta,
       spatial_scale = spatial_scale, step_scale = step_cells * spatial_scale)
}

# normalize compactness input for slic/metrics workflows
.sc_util_prep_compactness = function(compactness) {
  if (is.numeric(compactness) && length(compactness) == 1 && !is.na(compactness) && compactness > 0) {
    return(list(value = compactness, adaptive = FALSE, compactness_method = "constant"))
  }

  if (inherits(compactness, "sc_adaptive")) {
    return(list(value = 0, adaptive = TRUE, compactness_method = compactness$method))
  }
  stop("The 'compactness' argument must be a single positive number or use_adaptive()", call. = FALSE)
}
