# run the algorithm on the area defined by 'ext'
run_slic_chunks = function(ext, x, step, compactness, dist_name,
                           dist_fun, avg_fun_fun, avg_fun_name, clean,
                           iter, minarea, transform, input_centers, verbose,
                           diagnostics = FALSE, iter_diagnostics = FALSE){
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
                  iter = iter, minarea = minarea, input_centers = input_centers, verbose = verbose,
                  diagnostics = diagnostics, iter_diagnostics = iter_diagnostics)
  # returns the initial centers if iter = 0
  if (iter == 0){
    slic_sf = data.frame(stats::na.omit(slic[[2]]))
    slic_sf[["X1"]] = as.vector(terra::ext(x))[[1]] + (slic_sf[["X1"]] * terra::res(x)[[1]]) + (terra::res(x)[[1]]/2)
    slic_sf[["X2"]] = as.vector(terra::ext(x))[[4]] - (slic_sf[["X2"]] * terra::res(x)[[2]]) - (terra::res(x)[[1]]/2)
    slic_sf = sf::st_as_sf(slic_sf, coords = c("X1", "X2"))
    sf::st_crs(slic_sf) = terra::crs(x)
    return(slic_sf)
  }
  # transforms the output back to RGB color space if transform = "to_LAB"
  if (!is.null(transform)){
    if (transform == "to_LAB"){
      slic[[3]] = grDevices::convertColor(slic[[3]], from = "Lab", to = "sRGB") * 255
    }
  }
  # prepares the output: an sf object with supercells' ids, coordinates, and average values
  if (nrow(slic[[2]]) == 0 || all(slic[[2]] == 0)) stop("I cannot return supercells. This may be due to a large number of missing values in the 'x' object. Try to either trim your data to the non-NA area (e.g., with 'terra::trim()') or increase the number of expected supercells.", call. = FALSE)
  slic_sf = terra::rast(slic[[1]])
  terra::NAflag(slic_sf) = -1
  terra::crs(slic_sf) = terra::crs(x)
  terra::ext(slic_sf) = terra::ext(x)
  slic_sf = sf::st_as_sf(terra::as.polygons(slic_sf, dissolve = TRUE))
  if (nrow(slic_sf) > 0){
    empty_centers = slic[[2]][,1] != 0 | slic[[2]][,2] != 0
    slic_sf = cbind(slic_sf, stats::na.omit(slic[[2]][empty_centers, ]))
    names(slic_sf) = c("supercells", "x", "y", "geometry")
    slic_sf[["supercells"]] = slic_sf[["supercells"]] + 1
    slic_sf[["x"]] = as.vector(terra::ext(x))[[1]] + (slic_sf[["x"]] * terra::res(x)[[1]]) + (terra::res(x)[[1]]/2)
    slic_sf[["y"]] = as.vector(terra::ext(x))[[4]] - (slic_sf[["y"]] * terra::res(x)[[2]]) - (terra::res(x)[[1]]/2)
    colnames(slic[[3]]) = names(x)
    slic_sf = cbind(slic_sf, stats::na.omit(slic[[3]][empty_centers, , drop = FALSE]))
    slic_sf = suppressWarnings(sf::st_collection_extract(slic_sf, "POLYGON"))
    # slic_sf = sf::st_cast(slic_sf, "MULTIPOLYGON")
    if (diagnostics && !is.null(slic[[4]])) {
      attr(slic_sf, "diagnostics") = slic[[4]]
    }
    if (iter_diagnostics && !is.null(slic[[5]])) {
      attr(slic_sf, "iter_diagnostics") = slic[[5]]
    }
    return(slic_sf)
  }
}
