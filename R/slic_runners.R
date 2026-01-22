# Run SLIC on the raster chunk defined by 'ext'
run_slic_chunk_raster = function(ext, x, step, compactness, dist_name,
                                 dist_fun, avg_fun_fun, avg_fun_name, clean,
                                 iter, minarea, transform, input_centers, verbose,
                                 iter_diagnostics = FALSE){
  centers = TRUE
  if (is.character(x)){
    x = terra::rast(x)
  }

  # crops the input to the chunk extent
  x = x[ext[1]:ext[2], ext[3]:ext[4], drop = FALSE]

  # gets the number of rows and columns, and the values of the input
  mat = dim(x)[1:2]; mode(mat) = "integer"
  vals = terra::values(x, mat = TRUE, na.rm = FALSE)
  mode(vals) = "double"

  # transforms the input to LAB color space if transform = "to_LAB"
  if (!is.null(transform)){
    if (transform == "to_LAB"){
      if (ncol(vals) > 3) {
        vals = vals[, 1:3]
        warning("The provided raster has more than three layers: only the first three were kept for calculations", call. = FALSE)
      }
      vals = vals / 255
      vals = grDevices::convertColor(vals, from = "sRGB", to = "Lab")
    }
  }

  # runs the algorithm
  slic = run_slic(mat, vals = vals, step = step, compactness = compactness, clean = clean,
                  centers = centers, dist_name = dist_name, dist_fun = dist_fun,
                  avg_fun_fun = avg_fun_fun, avg_fun_name = avg_fun_name,
                  iter = iter, minarea = minarea, input_centers = input_centers,
                  verbose = verbose, iter_diagnostics = iter_diagnostics)

  # transforms the output back to RGB color space if transform = "to_LAB"
  if (!is.null(transform)){
    if (transform == "to_LAB"){
      slic[[3]] = grDevices::convertColor(slic[[3]], from = "Lab", to = "sRGB") * 255
    }
  }

  if (iter == 0){
    # returns the initial centers if iter = 0
    slic_sf = data.frame(stats::na.omit(slic[[2]]))
    slic_sf[["X1"]] = as.vector(terra::ext(x))[[1]] + (slic_sf[["X1"]] * terra::res(x)[[1]]) + (terra::res(x)[[1]]/2)
    slic_sf[["X2"]] = as.vector(terra::ext(x))[[4]] - (slic_sf[["X2"]] * terra::res(x)[[2]]) - (terra::res(x)[[2]]/2)
    slic_sf = sf::st_as_sf(slic_sf, coords = c("X1", "X2"))
    sf::st_crs(slic_sf) = terra::crs(x)
    result = list(iter0 = TRUE, raster = NULL, centers_sf = slic_sf,
                centers = slic[[2]], centers_vals = slic[[3]],
                iter_diagnostics = slic[[4]], names_x = names(x))
  } else {
    # prepares the output: an sf object with supercells' ids, coordinates, and average values
    if (nrow(slic[[2]]) == 0 || all(slic[[2]] == 0, na.rm = TRUE)) stop("I cannot return supercells. This may be due to a large number of missing values in the 'x' object. Try to either trim your data to the non-NA area (e.g., with 'terra::trim()') or increase the number of expected supercells.", call. = FALSE)
    slic[[1]] = slic[[1]] + 1
    slic_rast = terra::rast(slic[[1]])
    terra::NAflag(slic_rast) = 0
    terra::crs(slic_rast) = terra::crs(x)
    terra::ext(slic_rast) = terra::ext(x)
    result = list(iter0 = FALSE, raster = slic_rast, centers_sf = NULL,
                centers = slic[[2]], centers_vals = slic[[3]],
                iter_diagnostics = slic[[4]], names_x = names(x))
  }

  return(result)
}

# Run SLIC on a chunk and return point centers instead of polygons
run_slic_chunk_points = function(ext, x, step, compactness, dist_name,
                                 dist_fun, avg_fun_fun, avg_fun_name, clean,
                                 iter, minarea, transform, input_centers, verbose,
                                 iter_diagnostics = FALSE, metadata = TRUE){
  res = run_slic_chunk_raster(ext, x, step, compactness, dist_name,
                             dist_fun, avg_fun_fun, avg_fun_name, clean,
                             iter, minarea, transform, input_centers, verbose,
                             iter_diagnostics)
  if (res$iter0) {
    raster_ref = if (is.character(x)) terra::rast(x) else x
    raster_ref = raster_ref[ext[1]:ext[2], ext[3]:ext[4], drop = FALSE]
    points_sf = slic_centers_points(res$centers, raster_ref, res$centers_vals, res$names_x)
  } else {
    points_sf = slic_centers_points(res$centers, res$raster, res$centers_vals, res$names_x)
  }
  if (!isTRUE(metadata) && "x" %in% names(points_sf)) {
    points_sf = points_sf[, setdiff(names(points_sf), c("x", "y")), drop = FALSE]
  }
  if (iter_diagnostics && !is.null(res$iter_diagnostics)) {
    attr(points_sf, "iter_diagnostics") = res$iter_diagnostics
  }
  return(points_sf)
}

# Run SLIC on a chunk and return polygon supercells with attributes
run_slic_chunks = function(ext, x, step, compactness, dist_name,
                           dist_fun, avg_fun_fun, avg_fun_name, clean,
                           iter, minarea, transform, input_centers, verbose,
                           iter_diagnostics = FALSE, metadata = TRUE){
  res = run_slic_chunk_raster(ext, x, step, compactness, dist_name,
                              dist_fun, avg_fun_fun, avg_fun_name, clean,
                              iter, minarea, transform, input_centers, verbose,
                              iter_diagnostics)
  if (res$iter0) {
    return(res$centers_sf)
  }
  slic_sf = sf::st_as_sf(terra::as.polygons(res$raster, dissolve = TRUE))
  if (nrow(slic_sf) > 0){
    names(slic_sf)[1] = "supercells"
    centers_sf = slic_centers_points(res$centers, res$raster, res$centers_vals, res$names_x)
    centers_df = sf::st_drop_geometry(centers_sf)
    centers_df[["supercells"]] = NULL
    slic_sf = cbind(slic_sf, centers_df)
    if (!isTRUE(metadata)) {
      slic_sf = slic_sf[, setdiff(names(slic_sf), c("x", "y")), drop = FALSE]
    }
    slic_sf = suppressWarnings(sf::st_collection_extract(slic_sf, "POLYGON"))
    if (iter_diagnostics && !is.null(res$iter_diagnostics)) {
      attr(slic_sf, "iter_diagnostics") = res$iter_diagnostics
    }
    return(slic_sf)
  }
}

# Convert SLIC centers to an sf point layer with optional values
# Used by: run_slic_chunk_points and run_slic_chunks
slic_centers_points = function(centers, raster, centers_vals = NULL, names_x = NULL) {
  centers_df = as.data.frame(centers)
  empty_centers = centers_df[, 1] != 0 | centers_df[, 2] != 0
  centers_df = centers_df[empty_centers, , drop = FALSE]

  if (nrow(centers_df) == 0) {
    empty_sf = sf::st_sf(supercells = integer(), x = double(), y = double(),
                         geometry = sf::st_sfc(crs = terra::crs(raster)))
    return(empty_sf)
  }

  names(centers_df) = c("x", "y")
  ext = terra::ext(raster)
  res = terra::res(raster)
  centers_df[["x"]] = as.vector(ext)[[1]] + (centers_df[["x"]] * res[[1]]) + (res[[1]] / 2)
  centers_df[["y"]] = as.vector(ext)[[4]] - (centers_df[["y"]] * res[[2]]) - (res[[2]] / 2)
  centers_df[["supercells"]] = which(empty_centers)

  if (!is.null(centers_vals)) {
    if (!is.null(names_x)) {
      colnames(centers_vals) = names_x
    }
    centers_df = cbind(centers_df, centers_vals[empty_centers, , drop = FALSE])
  }

  centers_sf = sf::st_as_sf(centers_df, coords = c("x", "y"), remove = FALSE)
  sf::st_crs(centers_sf) = terra::crs(raster)
  
  return(centers_sf)
}
