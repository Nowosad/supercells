# check if raster is in memory
in_memory = function(x){
  terra::sources(x) == ""
}

# normalize input to SpatRaster
.sc_prep_raster = function(x) {
  if (inherits(x, "SpatRaster")) {
    return(x)
  }
  if (inherits(x, "stars")) {
    return(terra::rast(x))
  }
  stop("The SpatRaster class is expected as an input", call. = FALSE)
}

.sc_prep_dist_fun = function(dist_fun) {
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

.sc_prep_avg_fun = function(avg_fun) {
  if (is.character(avg_fun)) {
    avg_fun_name = avg_fun
    avg_fun_fun = function() ""
  } else {
    avg_fun_name = ""
    avg_fun_fun = avg_fun
  }
  return(list(avg_fun_name = avg_fun_name, avg_fun_fun = avg_fun_fun))
}

# converts sf object ('y') to a matrix of coordinates based on a raster ('x') dimensions
centers_to_dims = function(x, y){
  y_coords = sf::st_coordinates(sf::st_geometry(y))
  y_col = terra::colFromX(x, y_coords[, 1])
  y_row = terra::rowFromY(x, y_coords[, 2])
  center_dims = cbind(y_col, y_row)
  storage.mode(center_dims) = "integer"
  unique(center_dims)
}

# creates a sequence of integers from 'from' to 'to' with a step 'by' (including the 'to' value)
seq_last = function(from, to, by){
  vec = do.call(what = seq.int, args = list(from, to, by))
  if (utils::tail(vec, 1) != to) {
    return(c(vec, to))
  } else {
    return(vec)
  }
}
