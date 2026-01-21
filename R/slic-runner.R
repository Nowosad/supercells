# run the algorithm on the area defined by 'ext'
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
  vals = as.matrix(terra::as.data.frame(x, cells = FALSE, na.rm = FALSE))
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
  # returns the initial centers if iter = 0
  if (iter == 0){
    slic_sf = data.frame(stats::na.omit(slic[[2]]))
    slic_sf[["X1"]] = as.vector(terra::ext(x))[[1]] + (slic_sf[["X1"]] * terra::res(x)[[1]]) + (terra::res(x)[[1]]/2)
    slic_sf[["X2"]] = as.vector(terra::ext(x))[[4]] - (slic_sf[["X2"]] * terra::res(x)[[2]]) - (terra::res(x)[[1]]/2)
    slic_sf = sf::st_as_sf(slic_sf, coords = c("X1", "X2"))
    sf::st_crs(slic_sf) = terra::crs(x)
    return(list(iter0 = TRUE, raster = NULL, centers_sf = slic_sf,
                centers = slic[[2]], centers_vals = slic[[3]],
                iter_diagnostics = slic[[4]], names_x = names(x)))
  }
  # transforms the output back to RGB color space if transform = "to_LAB"
  if (!is.null(transform)){
    if (transform == "to_LAB"){
      slic[[3]] = grDevices::convertColor(slic[[3]], from = "Lab", to = "sRGB") * 255
    }
  }
  # prepares the output: an sf object with supercells' ids, coordinates, and average values
  if (nrow(slic[[2]]) == 0 || all(slic[[2]] == 0)) stop("I cannot return supercells. This may be due to a large number of missing values in the 'x' object. Try to either trim your data to the non-NA area (e.g., with 'terra::trim()') or increase the number of expected supercells.", call. = FALSE)
  slic_rast = terra::rast(slic[[1]])
  terra::NAflag(slic_rast) = -1
  terra::crs(slic_rast) = terra::crs(x)
  terra::ext(slic_rast) = terra::ext(x)
  list(iter0 = FALSE, raster = slic_rast, centers_sf = NULL,
       centers = slic[[2]], centers_vals = slic[[3]],
       iter_diagnostics = slic[[4]], names_x = names(x))
}

run_slic_chunks = function(ext, x, step, compactness, dist_name,
                           dist_fun, avg_fun_fun, avg_fun_name, clean,
                           iter, minarea, transform, input_centers, verbose,
                           iter_diagnostics = FALSE){
  res = run_slic_chunk_raster(ext, x, step, compactness, dist_name,
                              dist_fun, avg_fun_fun, avg_fun_name, clean,
                              iter, minarea, transform, input_centers, verbose,
                              iter_diagnostics)
  if (res$iter0) {
    return(res$centers_sf)
  }
  slic_sf = sf::st_as_sf(terra::as.polygons(res$raster, dissolve = TRUE))
  if (nrow(slic_sf) > 0){
    empty_centers = res$centers[,1] != 0 | res$centers[,2] != 0
    slic_sf = cbind(slic_sf, stats::na.omit(res$centers[empty_centers, ]))
    names(slic_sf) = c("supercells", "x", "y", "geometry")
    slic_sf[["supercells"]] = slic_sf[["supercells"]] + 1
    slic_sf[["x"]] = as.vector(terra::ext(res$raster))[[1]] + (slic_sf[["x"]] * terra::res(res$raster)[[1]]) + (terra::res(res$raster)[[1]]/2)
    slic_sf[["y"]] = as.vector(terra::ext(res$raster))[[4]] - (slic_sf[["y"]] * terra::res(res$raster)[[2]]) - (terra::res(res$raster)[[1]]/2)
    colnames(res$centers_vals) = res$names_x
    slic_sf = cbind(slic_sf, stats::na.omit(res$centers_vals[empty_centers, , drop = FALSE]))
    slic_sf = suppressWarnings(sf::st_collection_extract(slic_sf, "POLYGON"))
    if (iter_diagnostics && !is.null(res$iter_diagnostics)) {
      attr(slic_sf, "iter_diagnostics") = res$iter_diagnostics
    }
    return(slic_sf)
  }
}
