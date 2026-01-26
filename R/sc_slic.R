#' Create supercells
#'
#' Creates supercells from single- or multi-band rasters using an extended SLIC algorithm.
#' The function supports either a target number of supercells (`k`) or a fixed grid
#' spacing (`step`), as well as optional custom centers and chunked processing.
#'
#' @details
#' Use [`sc_slic()`] for polygon outputs. For raster or point centers outputs, see
#' [`sc_slic_raster()`] and [`sc_slic_points()`].
#' Evaluation and diagnostic options:
#' \itemize{
#'   \item Iteration diagnostics: set `iter_diagnostics = TRUE` to attach an
#'   `iter_diagnostics` attribute (only available without chunking). Use
#'   [`sc_plot_iter_diagnostics()`] to visualize the convergence over iterations.
#'   \item Pixel diagnostics: [sc_metrics_pixels()] for per-pixel spatial, value,
#'   and combined distances.
#'   \item Cluster diagnostics: [sc_metrics_supercells()] for per-supercell summaries.
#'   \item Global diagnostics: [sc_metrics_global()] for a single-row summary.
#' }
#' @seealso [`sc_slic_raster()`], [`sc_slic_points()`], [`sc_plot_iter_diagnostics()`],
#'   [`sc_metrics_pixels()`], [`sc_metrics_supercells()`], [`sc_metrics_global()`]
#'
#' @param x An object of class SpatRaster (terra) or class stars (stars).
#' @param step The distance (number of cells) between initial centers (alternative is `k`).
#' @param compactness A compactness value.
#' @param dist_fun A distance function name or a custom function. Supported names:
#' "euclidean", "jsd", "dtw", "dtw2d", or any method from `philentropy::getDistMethods()`.
#' A custom function must accept two numeric vectors and return a single numeric value.
#' @param avg_fun An averaging function name or custom function used to summarize
#' values within each supercell. Supported names: "mean" and "median". A custom
#' function must accept a numeric vector and return a single numeric value.
#' @param clean Should connectivity of the supercells be enforced?
#' @param minarea Minimal size of a supercell (in cells).
#' @param iter Number of iterations.
#' @param k The number of supercells desired (alternative to `step`).
#' @param centers Optional sf object of custom centers. Requires `step`.
#' @param metadata Logical. Should metadata columns be kept?
#' @param chunks Chunking option. Use `FALSE` for no chunking, `TRUE` for
#' automatic chunking based on size, or a numeric value for a fixed chunk size
#' (in number of cells per side).
#' @param future Logical. Use future for parallelization?
#' @param verbose Verbosity level.
#' @param iter_diagnostics Logical. If `TRUE`, attaches iteration diagnostics as an
#' attribute (`iter_diagnostics`) on the output. Only available when chunks are not used.
#'
#' @return An sf object with the supercell polygons and summary statistics.
#' Information on `step` and `compactness` are attached to the result as attributes.
#'
#' @references Achanta, R., Shaji, A., Smith, K., Lucchi, A., Fua, P., & Süsstrunk, S. (2012). SLIC Superpixels Compared to State-of-the-Art Superpixel Methods. IEEE Transactions on Pattern Analysis and Machine Intelligence, 34(11), 2274–2282. https://doi.org/10.1109/tpami.2012.120
#' @references Nowosad, J., Stepinski, T. (2022). Extended SLIC superpixels algorithm for applications to non-imagery geospatial rasters. International Journal of Applied Earth Observation and Geoinformation, https://doi.org/10.1016/j.jag.2022.102935
#' @export
#'
#' @examples
#' library(supercells)
#' # One variable
#'
#' vol = terra::rast(system.file("raster/volcano.tif", package = "supercells"))
#' vol_slic1 = sc_slic(vol, step = 8, compactness = 1)
#' terra::plot(vol)
#' plot(sf::st_geometry(vol_slic1), add = TRUE, lwd = 0.2)
sc_slic = function(x, step = NULL, compactness, dist_fun = "euclidean",
                   avg_fun = "mean", clean = TRUE, minarea, iter = 10,
                   k = NULL, centers = NULL, metadata = FALSE, chunks = FALSE,
                   future = FALSE, iter_diagnostics = FALSE, verbose = 0) {

  if (iter == 0) {
    stop("iter = 0 returns centers only; polygon output is not available. Use sc_slic_points(iter = 0) to get initial centers.", call. = FALSE)
  }

  prep_args = .sc_slic_prep_args(x, step, compactness, dist_fun, avg_fun, clean, minarea, iter,
                                 k, centers, metadata, chunks, future, iter_diagnostics, verbose)

  segment = .sc_slic_segment(prep_args, .sc_run_full_polygons, .sc_run_chunk_polygons)

  iter_attr = .sc_slic_add_iter_attr(segment$chunks, prep_args$iter_diagnostics)

  result = .sc_slic_post(segment$chunks, prep_args, iter_attr)
  return(result)
}
