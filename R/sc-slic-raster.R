#' Create supercells as a raster
#'
#' Runs the SLIC workflow and returns a raster of supercell IDs.
#' For polygon outputs, use `sc_slic`; for point centers, use `sc_slic_points`.
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
  prep = .sc_slic_prep(x, k, step, centers, compactness, dist_fun, avg_fun,
                       minarea, chunks, iter, transform, metadata, iter_diagnostics)

  if (iter == 0) {
    stop("iter = 0 returns centers only; raster output is not available", call. = FALSE)
  }

  if (nrow(prep$chunk_ext) == 1) {
    res_list = list(.sc_slic_run_single_raster(prep, compactness, clean, iter, transform, verbose, future))
  } else {
    res_list = .sc_slic_run_chunks_raster(prep, compactness, clean, iter, transform, verbose, future)
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

.sc_slic_run_single_raster = function(prep, compactness, clean, iter, transform, verbose, future) {
  ext = prep$chunk_ext[1, ]
  run_slic_chunk_raster(ext, prep$x, step = prep$step, compactness = compactness,
                        dist_name = prep$funs$dist_name, dist_fun = prep$funs$dist_fun,
                        avg_fun_fun = prep$funs$avg_fun_fun, avg_fun_name = prep$funs$avg_fun_name,
                        clean = clean, iter = iter, minarea = prep$minarea, transform = transform,
                        input_centers = prep$input_centers, verbose = verbose,
                        iter_diagnostics = prep$iter_diagnostics)
}

.sc_slic_run_chunks_raster = function(prep, compactness, clean, iter, transform, verbose, future) {
  .sc_slic_apply_chunks(prep, run_slic_chunk_raster, compactness, clean, iter, transform, verbose, future)
}
