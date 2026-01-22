#' Create supercells as a raster
#'
#' Runs the SLIC workflow and returns a raster of supercell IDs
#' IDs are 1-based and are unique across chunks when chunking is used
#' For polygon outputs, use `sc_slic`; for point centers, use `sc_slic_points`
#'
#' @inheritParams sc_slic
#'
#' @return A SpatRaster with supercell IDs
#' @export
#' @examples
#' library(supercells)
#' vol = terra::rast(system.file("raster/volcano.tif", package = "supercells"))
#' vol_ids = sc_slic_raster(vol, step = 8, compactness = 1)
#' terra::plot(vol_ids)
sc_slic_raster = function(x, step = NULL, compactness, dist_fun = "euclidean",
                          avg_fun = "mean", clean = TRUE, minarea, iter = 10,
                          transform = NULL, k = NULL, centers = NULL,
                          metadata = FALSE, chunks = FALSE, future = FALSE, verbose = 0,
                          iter_diagnostics = FALSE) {

  # prep arguments
  prep_args = .sc_slic_prep_args(x, step, compactness, k, centers, dist_fun, avg_fun,
                            minarea, chunks, iter, transform, metadata, iter_diagnostics)

  # run the exSLIC algorithm
  if (nrow(prep_args$chunk_ext) == 1) {
    res_list = .sc_slic_run_single_raster(prep_args, compactness, clean, iter, transform, verbose, future)
  } else {
    res_list = .sc_slic_run_chunks_raster(prep_args, compactness, clean, iter, transform, verbose, future)
  }

  # if (!is.list(res_list) || (is.list(res_list) && !is.list(res_list[[1]]))) {
  #   res_list = list(res_list)
  # }

  # update the supercells IDS and merge the results ()
  rasters = offset_supercells_ids_raster(lapply(res_list, function(x) x$raster))
  if (length(rasters) == 1) { # merge needs at least two rasters
    sc_rast = rasters[[1]]
  } else {
    sc_rast = do.call(terra::merge, rasters)
  }

  return(sc_rast)
}

.sc_slic_run_single_raster = function(prep, compactness, clean, iter, transform, verbose, future) {
  ext = prep$chunk_ext[1, ]
  result = list(run_slic_chunk_raster(ext, prep$x, step = prep$step, compactness = compactness,
                                      dist_name = prep$funs$dist_name, dist_fun = prep$funs$dist_fun,
                                      avg_fun_fun = prep$funs$avg_fun_fun, avg_fun_name = prep$funs$avg_fun_name,
                                      clean = clean, iter = iter, minarea = prep$minarea, transform = transform,
                                      input_centers = prep$input_centers, verbose = verbose,
                                      iter_diagnostics = prep$iter_diagnostics))
  return(result)
}

.sc_slic_run_chunks_raster = function(prep, compactness, clean, iter, transform, verbose, future) {
  return(.sc_slic_apply_chunks(prep, run_slic_chunk_raster, compactness, clean, iter, transform, verbose, future))
}
