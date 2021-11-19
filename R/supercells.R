#' Creates supercells
#'
#' Creates supercells based on single- or multi-band spatial raster data. It uses a modified version of the SLIC Superpixel algorithm by Achanta et al. (2012), allowing specification of a distance function.
#'
#' @param x An object of class SpatRaster (terra)
#' @param k The number of supercells desired by the user (the output number can be slightly different!).
#' You can use either `k` or `step`.
#' @param compactness A compactness value. Larger values cause clusters to be more compact/even (square).
#' A compactness value depends on the range of input cell values and selected distance measure.
#' @param dist_fun A distance function. Currently implemented distance functions are "euclidean", "jsd", "dtw" (dynamic time warping), name of any distance function from the `philentropy` package (see [philentropy::getDistMethods()]), or any user defined function accepting two vectors and returning one value. Default: "euclidean"
#' @param avg_fun An averaging function - how the values of the supercells' centers are calculated?
#' It accepts any fitting R function (e.g., `base::mean()` or `stats::median()`) or one of internally implemented `"mean"` and `"median"`. Default: `"mean"`
#' @param clean Should connectivity of the supercells be enforced?
#' @param iter The number of iterations performed to create the output.
#' @param minarea Specifies the minimal size of a supercell (in cells). Only works when `clean = TRUE`.
#' By default, when `clean = TRUE`, average area (A) is calculated based on the total number of cells divided by a number of superpixels.
#' Next, the minimal size of a supercell equals to A/(2^2) (A is being right shifted)
#' @param step The distance (number of cells) between initial supercells' centers. You can use either `k` or `step`.
#' @param transform Transformation to be performed on the input. Currently implemented is "to_LAB" allowing to convert RGB raster to a raster in the LAB color space. By default, no transformation is performed.
#' @param chunks Should the input (`x`) be split into chunks before deriving supercells? Either `FALSE` (default), `TRUE` (only large input objects are split), or a numeric value (representing the side length of the chunk in the number of cells).
#' @param future TRUE/FALSE
#'
#' @return An sf object with several columns: (1) supercells - an id of each supercell, (2) y and x coordinates, (3) one or more columns with average values of given variables in each supercell
#'
#' @references Achanta, R., Shaji, A., Smith, K., Lucchi, A., Fua, P., & Süsstrunk, S. (2012). SLIC Superpixels Compared to State-of-the-Art Superpixel Methods. IEEE Transactions on Pattern Analysis and Machine Intelligence, 34(11), 2274–2282. https://doi.org/10.1109/tpami.2012.120
#' @references Nowosad, J. Motif: an open-source R tool for pattern-based spatial analysis. Landscape Ecol (2021). https://doi.org/10.1007/s10980-020-01135-0
#' @export
#'
#' @examples
#' library(supercells)
#' library(terra)
#' library(sf)
#' # One variable
#'
#' vol = rast(system.file("raster/volcano.tif", package = "supercells"))
#' vol_slic1 = supercells(vol, k = 50, compactness = 1)
#' plot(vol)
#' plot(st_geometry(vol_slic1), add = TRUE, lwd = 0.2)
#'
#' # RGB variables
#'
#' ortho = rast(system.file("raster/ortho.tif", package = "supercells"))
#' ortho_slic1 = supercells(ortho, k = 1000, compactness = 10, transform = "to_LAB")
#' plot(ortho)
#' plot(st_geometry(ortho_slic1), add = TRUE)
#'
#' ### RGB variables - colored output
#'
#' rgb_to_hex = function(x){
#'   apply(t(x), 2, function(x) rgb(x[1], x[2], x[3], maxColorValue = 255))
#' }
#' avg_colors = rgb_to_hex(st_drop_geometry(ortho_slic1[4:6]))
#'
#' plot(ortho)
#' plot(st_geometry(ortho_slic1), add = TRUE, col = avg_colors)
supercells = function(x, k, compactness, dist_fun = "euclidean", avg_fun = "mean", clean = TRUE,
                      iter = 10, transform = NULL, step, minarea, chunks = FALSE, future = FALSE){
  if (!inherits(x, "SpatRaster")){
    stop("The SpatRaster class is expected as an input", call. = FALSE)
  }
  mat = dim(x)[1:2]
  mode(mat) = "integer"
  new_centers = matrix(c(0L, 0L), ncol = 2)
  if (!missing(k) && inherits(k, "sf")){
    new_centers = centers_to_dims(x, k)
  } else if (!missing(step) && !missing(k)){
    stop("You can specify either k or step, not both", call. = FALSE)
  } else if (missing(step) && missing(k)){
    stop("You need to specify either k or step", call. = FALSE)
  } else if (missing(step)){
    superpixelsize = round((mat[1] * mat[2]) / k + 0.5)
    step = round(sqrt(superpixelsize) + 0.5)
  }
  if (is.character(avg_fun)){
    avg_fun_name = avg_fun; avg_fun_fun = function() ""
  } else {
    avg_fun_name = ""; avg_fun_fun = avg_fun
  }
  if (is.character(dist_fun)){
    if (!(dist_fun %in% c("euclidean", "jsd", "dtw", philentropy::getDistMethods()))){
      stop("The provided distance function ('dist_fun') does not exist!", call. = FALSE)
    }
    dist_type = dist_fun; dist_fun = function() ""
  } else {
    dist_type = ""
  }
  if (missing(minarea)){
    minarea = 0
  }
  # split
  chunk_ext = prep_chunks_ext(dim(x), limit = chunks)
  if (!in_memory(x)){
    x = terra::sources(x)[["source"]][[1]]
  }
  if (future){
    oopts = options(future.globals.maxSize = +Inf)
    on.exit(options(oopts))
    slic_sf = future.apply::future_apply(chunk_ext, MARGIN = 1, run_slic_chunks, x = x,
                                         step = step, compactness = compactness, dist_type = dist_type,
                                         dist_fun = dist_fun, avg_fun_fun = avg_fun_fun, avg_fun_name = avg_fun_name,
                                         clean = clean, iter = iter, minarea = minarea, transform = transform, new_centers = new_centers, future.seed = TRUE)
  } else{
    slic_sf = apply(chunk_ext, MARGIN = 1, run_slic_chunks, x = x,
                                         step = step, compactness = compactness, dist_type = dist_type,
                                         dist_fun = dist_fun, avg_fun_fun = avg_fun_fun, avg_fun_name = avg_fun_name,
                                         clean = clean, iter = iter, minarea = minarea, transform = transform, new_centers = new_centers)
  }


  # combine
  slic_sf = update_supercells_ids(slic_sf)
  return(slic_sf)
}

# ext = c(1, 10, 1, 10)
run_slic_chunks = function(ext, x, step, compactness, dist_type,
                           dist_fun, avg_fun_fun, avg_fun_name, clean,
                           iter, minarea, transform, new_centers){
  centers = TRUE
  if (is.character(x)){
    x = terra::rast(x)
  }
  x = x[ext[1]:ext[2], ext[3]:ext[4], drop = FALSE]
  ext_x = terra::ext(x)
  mat = dim(x)[1:2]
  mode(mat) = "integer"
  vals = as.matrix(terra::as.data.frame(x, cell = FALSE, na.rm = FALSE))

  if (!is.null(transform)){
    if (transform == "to_LAB"){
      vals = vals / 255
      vals = grDevices::convertColor(vals, from = "sRGB", to = "Lab")
    }
  }
  slic = run_slic(mat, vals = vals, step = step, nc = compactness, con = clean,
                  centers = centers, type = dist_type, type_fun = dist_fun,
                  avg_fun_fun = avg_fun_fun, avg_fun_name = avg_fun_name,
                  iter = iter, lims = minarea, input_centers = new_centers)
  if (iter == 0){
    slic_sf = data.frame(stats::na.omit(slic[[2]]))
    slic_sf[["X1"]] = as.vector(terra::ext(x))[[1]] + (slic_sf[["X1"]] * terra::res(x)[[1]]) + (terra::res(x)[[1]]/2)
    slic_sf[["X2"]] = as.vector(terra::ext(x))[[4]] - (slic_sf[["X2"]] * terra::res(x)[[2]]) - (terra::res(x)[[1]]/2)
    slic_sf = sf::st_as_sf(slic_sf, coords = c("X1", "X2"))
    sf::st_crs(slic_sf) = terra::crs(x)
    return(slic_sf)
  }
  if (!is.null(transform)){
    if (transform == "to_LAB"){
      slic[[3]] = grDevices::convertColor(slic[[3]], from = "Lab", to = "sRGB") * 255
    }
  }
  slic_sf = terra::rast(slic[[1]])
  terra::NAflag(slic_sf) = -1
  terra::crs(slic_sf) = terra::crs(x)
  terra::ext(slic_sf) = ext_x
  slic_sf = sf::st_as_sf(terra::as.polygons(slic_sf, dissolve = TRUE))
  if (nrow(slic_sf) > 0){
    empty_centers = slic[[2]][,1] != 0 | slic[[2]][,2] != 0
    slic_sf = cbind(slic_sf, stats::na.omit(slic[[2]][empty_centers, ]))
    names(slic_sf) = c("supercells", "x", "y", "geometry")
    slic_sf[["supercells"]] = slic_sf[["supercells"]] + 1
    slic_sf[["x"]] = as.vector(ext_x)[[1]] + (slic_sf[["x"]] * terra::res(x)[[1]]) + (terra::res(x)[[1]]/2)
    slic_sf[["y"]] = as.vector(ext_x)[[4]] - (slic_sf[["y"]] * terra::res(x)[[2]]) - (terra::res(x)[[1]]/2)
    colnames(slic[[3]]) = names(x)
    slic_sf = cbind(slic_sf, stats::na.omit(slic[[3]][empty_centers, , drop = FALSE]))
    slic_sf = suppressWarnings(sf::st_collection_extract(slic_sf, "POLYGON"))
    slic_sf = sf::st_cast(slic_sf, "MULTIPOLYGON")
    return(slic_sf)
  }
}

update_supercells_ids = function(x){
  x = x[lapply(x, length) > 0]
  no_updates = length(x) - 1
  for (i in seq_len(no_updates)){
    prev_max = max(x[[i]][["supercells"]])
    x[[i + 1]][["supercells"]] = x[[i + 1]][["supercells"]] + prev_max
  }
  x = do.call(rbind, x)
  return(x)
}

pred_mem_usage = function(dim_x){
  mem_bytes = dim_x[1] * dim_x[2] * dim_x[3] * 8 #in bytes
  mem_gb = mem_bytes / (1024 * 1024 * 1024)
  mem_gb
}

optimize_chunk_size = function(dim_x, limit, by = 500){
  min_diff_memory = function(a, dim_x, limit){
    abs((dim_x[3] * a^2 * 8 / (1024 * 1024 * 1024)) - limit)
  }
  opti = stats::optimize(min_diff_memory,
                  interval = c(seq(100, max(dim_x[1:2]), by = by), max(dim_x[1:2])),
                  dim_x, limit)
  return(opti$minimum)
}

prep_chunks_ext = function(dim_x, limit){
  if (is.numeric(limit)){
    wsize = limit
    limit = 0
  } else if (!limit){
    limit = Inf
  } else {
    limit = 1 #hardcoded limit
    wsize = optimize_chunk_size(dim_x, limit, by = 500)
  }
  if (pred_mem_usage(dim_x) > limit){
    dims1 = ceiling(seq.int(0, to = dim_x[1],
                            length.out = as.integer((dim_x[1] - 1) / wsize + 1) + 1))
    dims2 = ceiling(seq.int(0, to = dim_x[2],
                            length.out = as.integer((dim_x[2] - 1) / wsize + 1) + 1))

    row_dims = seq_along(dims1)[-length(dims1)]
    col_dims = seq_along(dims2)[-length(dims2)]
    n_chunks = max(row_dims) * max(col_dims)
    row_cols_chunks = cbind(min_row = integer(length = n_chunks),
                            max_row = integer(length = n_chunks),
                            min_col = integer(length = n_chunks),
                            max_col = integer(length = n_chunks))
    l = 0
    for (i in row_dims){
      for (j in col_dims){
        l = l + 1
        row_cols_chunks[l, 1] = dims1[i] + 1
        row_cols_chunks[l, 2] = dims1[i + 1]
        row_cols_chunks[l, 3] = dims2[j] + 1
        row_cols_chunks[l, 4] = dims2[j + 1]
      }
    }
  } else{
    row_cols_chunks = cbind(min_row = 1,
                            max_row = dim_x[1],
                            min_col = 1,
                            max_col = dim_x[2])
  }
  return(row_cols_chunks)
}

in_memory = function(x){
  terra::sources(x)[["source"]] == ""
}

centers_to_dims = function(x, y){
  y_coords = sf::st_coordinates(sf::st_geometry(y))
  y_col = terra::colFromX(x, y_coords[, 1])
  y_row = terra::rowFromY(x, y_coords[, 2])
  center_dims = cbind(y_col, y_row)
  storage.mode(center_dims) = "integer"
  unique(center_dims)
}
