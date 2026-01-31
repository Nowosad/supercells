# SLIC runners and chunk helpers (internal)

# run slic on a chunk and return centers only (iter = 0)
.sc_run_chunk_centers = function(ext, x, step, compactness, dist_name,
                                 adaptive_compactness, dist_fun, avg_fun_fun, avg_fun_name, clean,
                                 minarea, input_centers,
                                 iter_diagnostics = FALSE, verbose = 0) {
  centers = TRUE
  if (is.character(x)) {
    x = terra::rast(x)
  }
  x = x[ext[1]:ext[2], ext[3]:ext[4], drop = FALSE]

  mat = dim(x)[1:2]; mode(mat) = "integer"
  vals = terra::values(x, mat = TRUE, na.rm = FALSE)
  mode(vals) = "double"

  slic = run_slic(mat = mat, vals = vals, step = step, compactness = compactness,
                  adaptive_compactness = adaptive_compactness, clean = clean, centers = centers,
                  dist_name = dist_name, dist_fun = dist_fun,
                  avg_fun_fun = avg_fun_fun, avg_fun_name = avg_fun_name,
                  iter = 0, minarea = minarea, input_centers = input_centers,
                  iter_diagnostics = iter_diagnostics, verbose = verbose)

  raster_ref = x
  list(centers = slic[[2]], centers_vals = slic[[3]],
       iter_diagnostics = slic[[4]], names_x = names(x),
       raster_ref = raster_ref)
}

# run slic on the raster chunk defined by ext
.sc_run_chunk_raster = function(ext, x, step, compactness, dist_name,
                                 adaptive_compactness, dist_fun, avg_fun_fun, avg_fun_name, clean,
                                 iter, minarea, input_centers,
                                 iter_diagnostics = FALSE, metadata = TRUE, verbose = 0){
  if (iter == 0) {
    stop("iter = 0 returns centers only; raster output is not available. Use sc_slic_points with iter = 0 to get initial centers.", call. = FALSE)
  }
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

  # runs the algorithm
  slic = run_slic(mat = mat, vals = vals, step = step, compactness = compactness,
                  adaptive_compactness = adaptive_compactness, clean = clean, centers = centers,
                  dist_name = dist_name, dist_fun = dist_fun,
                  avg_fun_fun = avg_fun_fun, avg_fun_name = avg_fun_name,
                  iter = iter, minarea = minarea, input_centers = input_centers,
                  iter_diagnostics = iter_diagnostics, verbose = verbose)

  # prepares the output: a raster of supercell ids plus centers and values
  if (nrow(slic[[2]]) == 0 || all(slic[[2]] == 0, na.rm = TRUE)) stop("I cannot return supercells. This may be due to a large number of missing values in the 'x' object. Try to either trim your data to the non-NA area (e.g., with 'terra::trim()') or increase the number of expected supercells.", call. = FALSE)
  slic[[1]] = slic[[1]] + 1
  slic_rast = x[[1]]
  terra::values(slic_rast) = slic[[1]]
  terra::NAflag(slic_rast) = 0
  result = list(iter0 = FALSE, raster = slic_rast, centers_sf = NULL,
                centers = slic[[2]], centers_vals = slic[[3]],
                iter_diagnostics = slic[[4]], names_x = names(x))

  return(result)
}

# convert slic centers to an sf point layer with optional values
# used by: .sc_run_chunk_points and .sc_run_chunk_polygons
.sc_run_centers_points = function(centers, raster, centers_vals = NULL, names_x = NULL) {
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
  if ("supercells" %in% names(centers_df)) {
    centers_df = centers_df[, c("supercells", setdiff(names(centers_df), "supercells")), drop = FALSE]
  }

  centers_sf = sf::st_as_sf(centers_df, coords = c("x", "y"), remove = FALSE)
  sf::st_crs(centers_sf) = terra::crs(raster)
  
  return(centers_sf)
}

# run slic on a chunk and return point centers instead of polygons
.sc_run_chunk_points = function(ext, x, step, compactness, dist_name,
                                 adaptive_compactness, dist_fun, avg_fun_fun, avg_fun_name, clean,
                                 iter, minarea, input_centers,
                                 iter_diagnostics = FALSE, metadata = TRUE, verbose = 0){
  if (iter == 0) {
    res = .sc_run_chunk_centers(ext, x, step, compactness, dist_name,
                               adaptive_compactness, dist_fun, avg_fun_fun, avg_fun_name, clean,
                               minarea, input_centers,
                               iter_diagnostics = iter_diagnostics, verbose = verbose)
    points_sf = .sc_run_centers_points(res$centers, res$raster_ref, res$centers_vals, res$names_x)
    # points_sf = stats::na.omit(points_sf)
    if (!isTRUE(metadata) && "x" %in% names(points_sf)) {
      points_sf = points_sf[, setdiff(names(points_sf), c("x", "y")), drop = FALSE]
    }
    if (iter_diagnostics && !is.null(res$iter_diagnostics)) {
      attr(points_sf, "iter_diagnostics") = res$iter_diagnostics
    }
    return(points_sf)
  }
  res = .sc_run_chunk_raster(ext, x, step, compactness, dist_name,
                             adaptive_compactness, dist_fun, avg_fun_fun, avg_fun_name, clean,
                             iter, minarea, input_centers,
                             iter_diagnostics = iter_diagnostics, verbose = verbose)
  points_sf = .sc_run_centers_points(res$centers, res$raster, res$centers_vals, res$names_x)
  ids = unique(terra::values(res$raster, mat = FALSE))
  ids = ids[!is.na(ids)]
  points_sf = points_sf[points_sf[["supercells"]] %in% ids, , drop = FALSE]
  points_sf = points_sf[match(ids, points_sf[["supercells"]]), , drop = FALSE]
  points_sf = stats::na.omit(points_sf)
  if (!isTRUE(metadata) && "x" %in% names(points_sf)) {
    points_sf = points_sf[, setdiff(names(points_sf), c("x", "y")), drop = FALSE]
  }
  if (iter_diagnostics && !is.null(res$iter_diagnostics)) {
    attr(points_sf, "iter_diagnostics") = res$iter_diagnostics
  }
  return(points_sf)
}

# run slic on a chunk and return polygon supercells with attributes
.sc_run_chunk_polygons = function(ext, x, step, compactness, dist_name,
                           adaptive_compactness, dist_fun, avg_fun_fun, avg_fun_name, clean,
                           iter, minarea, input_centers,
                           iter_diagnostics = FALSE, metadata = TRUE, verbose = 0){
  if (iter == 0) {
    stop("iter = 0 returns centers only; polygon output is not available. Use sc_slic_points with iter = 0 to get initial centers.", call. = FALSE)
  }
  res = .sc_run_chunk_raster(ext, x, step, compactness, dist_name,
                              adaptive_compactness, dist_fun, avg_fun_fun, avg_fun_name, clean,
                              iter, minarea, input_centers,
                              iter_diagnostics = iter_diagnostics, verbose = verbose)
  slic_sf = sf::st_as_sf(terra::as.polygons(res$raster, dissolve = TRUE))
  if (nrow(slic_sf) > 0){
    names(slic_sf)[1] = "supercells"
    centers_sf = .sc_run_centers_points(res$centers, res$raster, res$centers_vals, res$names_x)
    centers_df = sf::st_drop_geometry(centers_sf)
    if ("supercells" %in% names(centers_df)) {
      centers_df = centers_df[centers_df[["supercells"]] %in% slic_sf[["supercells"]], , drop = FALSE]
      centers_df = centers_df[match(slic_sf[["supercells"]], centers_df[["supercells"]]), , drop = FALSE]
      centers_df[["supercells"]] = NULL
    }
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

# run slic on a full raster using a chunk runner
.sc_run_full = function(x, runner, ...) {
  if (is.character(x)) {
    x = terra::rast(x)
  }
  ext = c(1L, dim(x)[1], 1L, dim(x)[2])
  runner(ext, x, ...)
}

# run slic on a full raster and return point centers
.sc_run_full_points = function(x, step, compactness, dist_name, dist_fun,
                           adaptive_compactness, avg_fun_fun, avg_fun_name, clean, iter, minarea,
                           input_centers, iter_diagnostics = FALSE,
                           metadata = TRUE, verbose = 0) {
  .sc_run_full(x, .sc_run_chunk_points, step, compactness, dist_name,
                    adaptive_compactness, dist_fun, avg_fun_fun, avg_fun_name, clean, iter, minarea,
                    input_centers, iter_diagnostics, metadata, verbose)
}

# run slic on a full raster and return raster ids
.sc_run_full_raster = function(x, step, compactness, dist_name, dist_fun,
                           adaptive_compactness, avg_fun_fun, avg_fun_name, clean, iter, minarea,
                           input_centers, iter_diagnostics = FALSE,
                           metadata = TRUE, verbose = 0) {
  .sc_run_full(x, .sc_run_chunk_raster, step, compactness, dist_name,
                    adaptive_compactness, dist_fun, avg_fun_fun, avg_fun_name, clean, iter, minarea,
                    input_centers, iter_diagnostics, metadata, verbose)
}

# run slic on a full raster and return polygon supercells
.sc_run_full_polygons = function(x, step, compactness, dist_name, dist_fun,
                             adaptive_compactness, avg_fun_fun, avg_fun_name, clean, iter, minarea,
                             input_centers, iter_diagnostics = FALSE,
                             metadata = TRUE, verbose = 0) {
  .sc_run_full(x, .sc_run_chunk_polygons, step, compactness, dist_name,
                    adaptive_compactness, dist_fun, avg_fun_fun, avg_fun_name, clean, iter, minarea,
                    input_centers, iter_diagnostics, metadata, verbose)
}
