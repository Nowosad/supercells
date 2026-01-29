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
                          step_unit = "cells", k = NULL, centers = NULL, metadata = FALSE,
                          chunks = FALSE,
                          iter_diagnostics = FALSE, verbose = 0) {

  if (iter == 0) {
    stop("iter = 0 returns centers only; raster output is not available. Use sc_slic_points(iter = 0) to get initial centers.", call. = FALSE)
  }
  # prep arguments
  prep_args = .sc_slic_prep_args(x, step, step_unit, compactness, dist_fun, avg_fun, clean, minarea, iter,
                            k, centers, metadata, chunks, iter_diagnostics, verbose)

  # segment once (single) or per chunk (chunked), returning a list of chunk results
  if (nrow(prep_args$chunk_ext) > 1) {
    n_chunks = nrow(prep_args$chunk_ext)
    expected_ids = .sc_chunk_expected_ids(prep_args$chunk_ext, prep_args$step)
    offsets = .sc_chunk_offsets(prep_args$chunk_ext, prep_args$step)
    max_expected = if (all(is.na(expected_ids))) NA_real_ else max(offsets + expected_ids, na.rm = TRUE)
    dtype = .sc_chunk_id_datatype(max_expected)
    chunk_files = character(n_chunks)
    on.exit(unlink(chunk_files), add = TRUE)
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
                                 prep_args$metadata, prep_args$verbose_cpp)
      r = res[["raster"]]
      if (!is.na(offsets[i]) && offsets[i] > 0) {
        r = r + offsets[i]
      }
      chunk_files[i] = tempfile(fileext = ".tif")
      wopt = list(gdal = c("TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256", "COMPRESS=NONE"))
      if (is.null(dtype)) {
        terra::writeRaster(r, filename = chunk_files[i], overwrite = TRUE, wopt = wopt)
      } else {
        terra::writeRaster(r, filename = chunk_files[i], overwrite = TRUE, wopt = wopt, datatype = dtype)
      }
    }
    out_file = tempfile(fileext = ".tif")
    if (is.numeric(prep_args$verbose) && prep_args$verbose > 0) {
      message("Merging chunk rasters...")
    }
    wopt = list(gdal = c("TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256", "COMPRESS=NONE"))
    if (!is.null(dtype)) {
      wopt[["datatype"]] = dtype
    }
    merged = terra::merge(terra::sprc(chunk_files), filename = out_file, overwrite = TRUE, wopt = wopt)
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
