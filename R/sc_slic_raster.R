#' Create supercells as a raster
#'
#' Runs the SLIC workflow and returns a raster of supercell IDs
#' IDs are 1-based and are unique across chunks when chunking is used
#' For polygon outputs, use [`sc_slic()`]; for point centers, use [`sc_slic_points()`]
#'
#' @inheritParams sc_slic
#' @param outcomes Character vector controlling which fields are returned.
#' Only `"supercells"` is supported in `sc_slic_raster()`.
#' @seealso [`sc_slic()`], [`sc_slic_points()`]
#'
#' @return A SpatRaster with supercell IDs.
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
                          k = NULL, centers = NULL,
                          outcomes = "supercells", chunks = FALSE,
                          verbose = 0) {

  if (iter == 0) {
    stop("iter = 0 returns centers only; raster output is not available. Use sc_slic_points(iter = 0) to get initial centers.", call. = FALSE)
  }
  if (!identical(outcomes, "supercells")) {
    stop("sc_slic_raster() supports only outcomes = 'supercells'", call. = FALSE)
  }
  # prep arguments
  prep_args = .sc_slic_prep_args(x, step, compactness, dist_fun, avg_fun, clean, minarea, iter,
                            k, centers, outcomes, chunks, verbose)

  # segment once (single) or per chunk (chunked), returning a list of chunk results
  if (nrow(prep_args$chunk_ext) > 1) {
    n_chunks = nrow(prep_args$chunk_ext)
    expected_ids = .sc_chunk_expected_max_ids(prep_args$chunk_ext, prep_args$step)
    max_expected = sum(expected_ids)
    dtype = .sc_chunk_id_datatype(max_expected)
    wopt = list(gdal = c("TILED=YES", "BLOCKXSIZE=256", "BLOCKYSIZE=256", "COMPRESS=NONE"), datatype = dtype)
    chunk_files = character(n_chunks)
    on.exit(unlink(chunk_files), add = TRUE)
    max_id = 0
    for (i in seq_len(n_chunks)) {
      if (is.numeric(prep_args$verbose) && prep_args$verbose > 0) {
        message(sprintf("Processing chunk %d/%d", i, n_chunks))
      }
      ext = prep_args$chunk_ext[i, ]
      res = .sc_run_chunk_raster(
        ext = ext,
        x = prep_args$x,
        step = prep_args$step,
        compactness = prep_args$compactness,
        dist_name = prep_args$funs$dist_name,
        adaptive_compactness = prep_args$adaptive_compactness,
        dist_fun = prep_args$funs$dist_fun,
        avg_fun_fun = prep_args$funs$avg_fun_fun,
        avg_fun_name = prep_args$funs$avg_fun_name,
        clean = prep_args$clean,
        iter = prep_args$iter,
        minarea = prep_args$minarea,
        input_centers = prep_args$input_centers,
        iter_diagnostics = FALSE,
        verbose = prep_args$verbose_cpp
      )
      r = res[["raster"]]
      if (max_id > 0) {
        r = r + max_id
      }
      n_centers = nrow(res[["centers"]])
      if (n_centers > 0) {
        max_id = max_id + n_centers
      }
      chunk_files[i] = tempfile(fileext = ".tif")
      terra::writeRaster(r, filename = chunk_files[i], overwrite = TRUE, wopt = wopt)
    }
    out_file = tempfile(fileext = ".tif")
    if (is.numeric(prep_args$verbose) && prep_args$verbose > 0) {
      message("Merging chunk rasters...")
    }
    result = terra::merge(terra::sprc(chunk_files), filename = out_file, overwrite = TRUE, wopt = wopt)
  } else {
    segment = .sc_slic_segment(prep_args, .sc_run_full_raster, .sc_run_chunk_raster)
    # single chunk: offset ids and return the raster
    chunk_rasters = .sc_chunk_offset_ids_raster_by_centers(segment$chunks)
    result = chunk_rasters[[1]]
  }
  names(result) = "supercells"
  return(result)
}
