#' Calculates initial maximum distances
#'
#' Calculates initial maximum distances from each supercell's center to the other cells. This includes maximum value distance, maximum spatial distance, and maximum total distance. This function can be useful to determine the initial parameters (especially the compactness value) for the supercells algorithm.
#'
#' @param x An object of class SpatRaster (terra) or class stars (stars)
#' @param k The number of supercells desired by the user (the output number can be slightly different!).
#' You can use either `k` or `step`.
#' It is also possible to provide a set of points (an `sf` object) as `k` together with the `step` value to create custom cluster centers.
#' @param compactness A compactness value. Larger values cause clusters to be more compact/even (square).
#' A compactness value depends on the range of input cell values and selected distance measure. By default, 1.
#' @param dist_fun A distance function. Currently implemented distance functions are `"euclidean"`, `"jsd"`, `"dtw"` (dynamic time warping), name of any distance function from the `philentropy` package (see [philentropy::getDistMethods()]; "log2" is used in this case), or any user defined function accepting two vectors and returning one value. Default: `"euclidean"`
#' @param avg_fun An averaging function - how the values of the supercells' centers are calculated?
#' It accepts any fitting R function (e.g., `base::mean()` or `stats::median()`) or one of internally implemented `"mean"` and `"median"`. Default: `"mean"`
#' @param step The distance (number of cells) between initial supercells' centers. You can use either `k` or `step`.
#' @param transform Transformation to be performed on the input. By default, no transformation is performed. Currently available transformation is "to_LAB": first, the conversion from RGB to the LAB color space is applied, then the supercells algorithm is run, and afterward, a reverse transformation is performed on the obtained results. (This argument is experimental and may be removed in the future).
#' @param metadata Logical. If `TRUE`, the output object will have metadata columns ("supercells", "x", "y"). If `FALSE`, the output object will not have metadata columns.
#' @param chunks Should the input (`x`) be split into chunks before deriving supercells? Either `FALSE` (default), `TRUE` (only large input objects are split), or a numeric value (representing the side length of the chunk in the number of cells).
#' @param future Should the future package be used for parallelization of the calculations? Default: `FALSE`. If `TRUE`, you also need to specify `future::plan()`.
#' @param verbose An integer specifying the level of text messages printed during calculations. 0 means no messages (default), 1 provides basic messages (e.g., calculation stage).
#'
#' @return An sf object with several columns: (1) supercells - an id of each supercell, (2) y and x coordinates, (3) one or more columns with average values of given variables in each supercell, (4) three columns with the maximum distances (`max_value_dist` -- maximum value distance, `max_spatial_dist` -- maximum spatial distance, `max_total_dist` -- maximum total distance) from a supercell's center to the other cells.
#'
#' @references Achanta, R., Shaji, A., Smith, K., Lucchi, A., Fua, P., & Süsstrunk, S. (2012). SLIC Superpixels Compared to State-of-the-Art Superpixel Methods. IEEE Transactions on Pattern Analysis and Machine Intelligence, 34(11), 2274–2282. https://doi.org/10.1109/tpami.2012.120
#' @references Nowosad, J. Motif: an open-source R tool for pattern-based spatial analysis. Landscape Ecol (2021). https://doi.org/10.1007/s10980-020-01135-0
#' @export
#'
#' @seealso [supercells()]
#'
#' @examples
#' library(supercells)
#' # One variable
#'
#' vol = terra::rast(system.file("raster/volcano.tif", package = "supercells"))
#' vol_dists = get_initial_distances(vol, k = 50, compactness = 1)
#' vol_dists
#' # plot(vol_dists)
#'
#' # RGB variables
#' # ortho = terra::rast(system.file("raster/ortho.tif", package = "supercells"))
#' # get_initial_distances(ortho, k = 1000, compactness = 1, metadata = FALSE)
get_initial_distances = function(x, k, compactness = 1, dist_fun = "euclidean", avg_fun = "mean", transform = NULL, step, metadata = TRUE, chunks = FALSE, future = FALSE, verbose = 0){
  result = supercells(x = x, k = k, compactness = c(compactness, 0),
                      dist_fun = dist_fun, avg_fun = avg_fun, clean = FALSE,
                      iter = 1, transform = transform, step = step,
                      metadata = metadata, chunks = chunks, future = future, verbose = verbose)
  return(result)
}

