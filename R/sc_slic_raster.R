#' Create supercells as a raster
#'
#' Runs the SLIC workflow and returns a raster of supercell IDs
#' IDs are 1-based and are unique across chunks when chunking is used
#' For polygon outputs, use [`sc_slic()`]; for point centers, use [`sc_slic_points()`]
#'
#' @inheritParams sc_slic
#' @seealso [`sc_slic()`]
#'
#' @return A SpatRaster with supercell IDs
#' 
#' @export
#' 
#' @examples
#' library(supercells)
#' vol = terra::rast(system.file("raster/volcano.tif", package = "supercells"))
#' vol_ids = sc_slic_raster(vol, step = 8, compactness = 1)
#' terra::plot(vol_ids)
sc_slic_raster = function(x, step = NULL, compactness, dist_fun = "euclidean",
                          avg_fun = "mean", clean = TRUE, minarea, iter = 10,
                          k = NULL, centers = NULL, metadata = FALSE,
                          chunks = FALSE, future = FALSE, iter_diagnostics = FALSE,
                          verbose = 0) {

  if (iter == 0) {
    stop("iter = 0 returns centers only; raster output is not available. Use sc_slic_points(iter = 0) to get initial centers.", call. = FALSE)
  }

  # prep arguments
  prep_args = .sc_slic_prep_args(x, step, compactness, dist_fun, avg_fun, clean, minarea, iter,
                            k, centers, metadata, chunks, future, iter_diagnostics, verbose)

  # segment once (single) or per chunk (chunked), returning a list of chunk results
  segment = .sc_slic_segment(prep_args, .sc_run_full_raster, .sc_run_chunk_raster)
  
  # merge/offset chunk rasters into a single ID raster
  chunk_rasters = lapply(segment$chunks, `[[`, "raster")
  chunk_rasters = .sc_chunk_offset_ids_raster(chunk_rasters)

  if (length(chunk_rasters) == 1) {
    result = chunk_rasters[[1]]
    names(result) = "supercells"
    return(result)
  } else {
    stitched = do.call(terra::merge, chunk_rasters)
    names(stitched) = "supercells"
    return(stitched)
  }
}
