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
                          chunks = FALSE,
                          iter_diagnostics = FALSE, verbose = 0) {

  if (iter == 0) {
    stop("iter = 0 returns centers only; raster output is not available. Use sc_slic_points(iter = 0) to get initial centers.", call. = FALSE)
  }

  # prep arguments
  prep_args = .sc_slic_prep_args(x, step, compactness, dist_fun, avg_fun, clean, minarea, iter,
                            k, centers, metadata, chunks, iter_diagnostics, verbose)

  # segment once (single) or per chunk (chunked), returning a list of chunk results
  if (nrow(prep_args$chunk_ext) > 1) {
    max_id = 0
    n_chunks = nrow(prep_args$chunk_ext)
    chunk_files = character(n_chunks)
    for (i in seq_len(n_chunks)) {
      if (is.numeric(prep_args$verbose) && prep_args$verbose > 0) {
        message(sprintf("Processing chunk %d/%d", i, n_chunks))
      }
      ext = prep_args$chunk_ext[i, ]
      res = .sc_run_chunk_raster(ext, prep_args$x, prep_args$step, prep_args$compactness,
                                 prep_args$funs$dist_name, prep_args$funs$dist_fun,
                                 prep_args$funs$avg_fun_fun, prep_args$funs$avg_fun_name,
                                 prep_args$clean, prep_args$iter, prep_args$minarea,
                                 prep_args$input_centers, prep_args$iter_diagnostics,
                                 prep_args$metadata, prep_args$verbose)
      r = res[["raster"]]
      if (max_id > 0) {
        r = r + max_id
      }
      curr_max = terra::global(r, "max", na.rm = TRUE)[1, 1]
      if (!is.na(curr_max)) {
        max_id = curr_max + 1
      }
      chunk_files[i] = tempfile(fileext = ".tif")
      terra::writeRaster(r, chunk_files[i], overwrite = TRUE)
    }
    out_file = tempfile(fileext = ".tif")
    merged = terra::merge(terra::sprc(chunk_files), filename = out_file, overwrite = TRUE)
    on.exit(unlink(chunk_files), add = TRUE)
    names(merged) = "supercells"
    return(merged)
  } else {
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
}
