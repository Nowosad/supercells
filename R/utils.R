# check if raster is in memory
in_memory = function(x){
  terra::sources(x) == ""
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
