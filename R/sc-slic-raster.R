#' Create supercells as a raster
#'
#' Runs [sc_slic()] and rasterizes the resulting supercells.
#'
#' @inheritParams sc_slic
#'
#' @return A SpatRaster with supercell IDs.
#' @export
sc_slic_raster = function(x, k = NULL, step = NULL, centers = NULL, compactness,
                          dist_fun = "euclidean", avg_fun = "mean", clean = TRUE,
                          iter = 10, transform = NULL, minarea, metadata = FALSE,
                          chunks = FALSE, future = FALSE, verbose = 0,
                          iter_diagnostics = FALSE) {
  x = .sc_create_input(x)
  .sc_create_validate(k, step, centers, compactness, chunks)

  if (iter == 0) {
    stop("iter = 0 returns centers only; raster output is not available", call. = FALSE)
  }
  step = .sc_create_step(x, k, step)
  input_centers = .sc_create_centers(x, centers)
  funs = .sc_create_funs(dist_fun, avg_fun)
  minarea = .sc_create_minarea(minarea, step)
  chunk_ext = prep_chunks_ext(dim(x), limit = chunks)
  if (iter_diagnostics && nrow(chunk_ext) > 1) {
    warning("Iteration diagnostics are only available when chunks = FALSE (single chunk). Iteration diagnostics were disabled.", call. = FALSE)
    iter_diagnostics = FALSE
  }

  if (future) {
    if (in_memory(x)) {
      names_x = names(x)
      x = terra::writeRaster(x, tempfile(fileext = ".tif"))
      names(x) = names_x
    }
    if (!in_memory(x)) {
      x = terra::sources(x)[[1]]
    }
    oopts = options(future.globals.maxSize = +Inf)
    on.exit(options(oopts))
    res_list = future.apply::future_apply(chunk_ext, MARGIN = 1, run_slic_chunk_raster, x = x,
                                          step = step, compactness = compactness, dist_name = funs$dist_name,
                                          dist_fun = funs$dist_fun, avg_fun_fun = funs$avg_fun_fun,
                                          avg_fun_name = funs$avg_fun_name, clean = clean, iter = iter,
                                          minarea = minarea, transform = transform, input_centers = input_centers,
                                          verbose = verbose, iter_diagnostics = iter_diagnostics,
                                          future.seed = TRUE)
  } else {
    res_list = apply(chunk_ext, MARGIN = 1, run_slic_chunk_raster, x = x,
                     step = step, compactness = compactness, dist_name = funs$dist_name,
                     dist_fun = funs$dist_fun, avg_fun_fun = funs$avg_fun_fun,
                     avg_fun_name = funs$avg_fun_name, clean = clean, iter = iter,
                     minarea = minarea, transform = transform, input_centers = input_centers,
                     verbose = verbose, iter_diagnostics = iter_diagnostics)
  }

  if (is.list(res_list) && length(res_list) > 0 && !is.null(res_list[[1]]$iter0) && res_list[[1]]$iter0) {
    stop("iter = 0 returns centers only; raster output is not available", call. = FALSE)
  }

  if (!is.list(res_list) || (is.list(res_list) && !is.list(res_list[[1]]))) {
    res_list = list(res_list)
  }

  rasters = lapply(res_list, function(x) x$raster)
  max_id = 0
  for (i in seq_along(rasters)) {
    r = rasters[[i]]
    if (max_id > 0) {
      r = terra::ifel(is.na(r), NA, r + max_id)
    }
    curr_max = terra::global(r, "max", na.rm = TRUE)[1, 1]
    if (!is.na(curr_max)) {
      max_id = max_id + curr_max
    }
    rasters[[i]] = r
  }
  sc_rast = Reduce(terra::merge, rasters)
  sc_rast = terra::ifel(is.na(sc_rast), NA, sc_rast + 1)

  sc_rast
}
