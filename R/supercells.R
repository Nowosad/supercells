#' Creates supercells
#'
#' Creates supercells by bringing the idea of superpixels into spatial analysis and grouping raster cells into coherent regions.
#' It uses a modified version of the SLIC algorithm by Achanta et al. (2012), allowing specification of a distance function.
#'
#' @param x An object of class SpatRaster (terra) or class stars (stars)
#' @param k The number of supercells desired by the user (the output number can be slightly different!).
#' You can use either `k` or `step`.
#' It is also possible to provide a set of points (an `sf` object) as `k` together with the `step` value to create custom cluster centers.
#' @param compactness A compactness value. Larger values cause supercells to be more compact/even (square).
#' Set `compactness = "auto"` to enable adaptive compactness (SLIC0).
#' A compactness value depends on the range of input cell values and selected distance measure.
#' @param dist_fun A distance function. Currently implemented distance functions are `"euclidean"`, `"jsd"`, `"dtw"` (dynamic time warping), name of any distance function from the `philentropy` package (see [philentropy::getDistMethods()]; "log2" is used in this case), or any user defined function accepting two vectors and returning one value. Default: `"euclidean"`
#' @param avg_fun An averaging function - how the values of the supercells' centers are calculated? The algorithm internally implements common functions `"mean"` and `"median"` (provided with quotation marks), but also accepts any fitting R function (e.g., `base::mean()` or `stats::median()`, provided as plain function name: `mean`). Default: `"mean"`. See details for more information.
#' @param clean Should connectivity of the supercells be enforced?
#' @param iter The number of iterations performed to create the output.
#' @param minarea Specifies the minimal size of a supercell (in cells). Only works when `clean = TRUE`.
#' By default, when `clean = TRUE`, average area (A) is calculated based on the total number of cells divided by a number of supercells
#' Next, the minimal size of a supercell equals to A/(2^2) (A is being right shifted)
#' @param step Initial center spacing. You can use either `k` or `step`.
#' Provide a plain numeric value for cell units, or use [in_meters()] for
#' map-distance steps in meters (automatically converted to cells using raster resolution).
#' @param transform Transformation to be performed on the input. By default, no transformation is performed. Currently available transformation is "to_LAB": first, the conversion from RGB to the LAB color space is applied, then the supercells algorithm is run, and afterward, a reverse transformation is performed on the obtained results. (This argument is experimental and may be removed in the future).
#' @param metadata Logical. Controls whether metadata columns
#' ("supercells", "x", "y") are included.
#' @param chunks Should the input (`x`) be processed with chunking before deriving supercells?
#' Use `FALSE` (default), `TRUE` (automatic chunking based on size), or a numeric value
#' (chunk side length in number of cells). When `TRUE`, the memory limit can be set
#' with `options(supercells.chunk_mem_gb = 4)`.
#' @param verbose An integer specifying the level of text messages printed during calculations. 0 means no messages (default), 1 provides basic messages (e.g., calculation stage).
#'
#' @details
#' If you want to use additional arguments for the averaging function (`avg_fun`), you can create a custom function. For example, if you want to calculate the mean by removing missing values, you can use the following code: `my_mean = function(x) mean(x, na.rm = TRUE)` and then provide `avg_fun = my_mean.`
#' For raster IDs or point centers outputs, see [sc_slic_raster()] and
#' [sc_slic_points()]. For evaluation and diagnostics, see
#' [sc_metrics_pixels()], [sc_metrics_supercells()], and [sc_metrics_global()].
#'
#' @return An sf object with several columns: (1) supercells - an id of each supercell, (2) y and x coordinates, (3) one or more columns with average values of given variables in each supercell.
#'
#' @references Achanta, R., Shaji, A., Smith, K., Lucchi, A., Fua, P., & Süsstrunk, S. (2012). SLIC Superpixels Compared to State-of-the-Art Superpixel Methods. IEEE Transactions on Pattern Analysis and Machine Intelligence, 34(11), 2274–2282. https://doi.org/10.1109/tpami.2012.120
#' @references Nowosad, J. Motif: an open-source R tool for pattern-based spatial analysis. Landscape Ecol (2021). https://doi.org/10.1007/s10980-020-01135-0
#' @seealso [`sc_slic()`]
#' @export
#'
#' @examples
#' library(supercells)
#' # One variable
#'
#' vol = terra::rast(system.file("raster/volcano.tif", package = "supercells"))
#' vol_slic1 = supercells(vol, k = 50, compactness = 1)
#' terra::plot(vol)
#' plot(sf::st_geometry(vol_slic1), add = TRUE, lwd = 0.2)
#'
#' # RGB variables
#' # ortho = terra::rast(system.file("raster/ortho.tif", package = "supercells"))
#' # ortho_slic1 = supercells(ortho, k = 1000, compactness = 10, transform = "to_LAB")
#' # terra::plot(ortho)
#' # plot(sf::st_geometry(ortho_slic1), add = TRUE)
#' #
#' # ### RGB variables - colored output
#' #
#' # rgb_to_hex = function(x){
#' #   apply(t(x), 2, function(x) rgb(x[1], x[2], x[3], maxColorValue = 255))
#' # }
#' # avg_colors = rgb_to_hex(sf::st_drop_geometry(ortho_slic1[4:6]))
#' #
#' # terra::plot(ortho)
#' # plot(sf::st_geometry(ortho_slic1), add = TRUE, col = avg_colors)
supercells = function(x, k, compactness, dist_fun = "euclidean", avg_fun = "mean", clean = TRUE,
                      iter = 10, transform = NULL, step, minarea, metadata = TRUE,
                      chunks = FALSE, verbose = 0){
  if (!is.logical(metadata) || length(metadata) != 1 || is.na(metadata)) {
    stop("The 'metadata' argument must be TRUE or FALSE", call. = FALSE)
  }
  outcomes = if (isTRUE(metadata)) {
    c("supercells", "coordinates", "values")
  } else {
    "values"
  }

  centers_arg = NULL
  if (!missing(k) && inherits(k, "sf")) {
    centers_arg = k
    k = NULL
  }

  x = .sc_util_prep_raster(x)

  trans = .supercells_transform_to_lab(x, transform)
  x = trans$x

  args = list(
    x,
    compactness = compactness,
    dist_fun = dist_fun,
    avg_fun = avg_fun,
    clean = clean,
    iter = iter,
    outcomes = outcomes,
    chunks = chunks,
    iter_diagnostics = FALSE,
    verbose = verbose
  )
  if (!missing(step)) {
    args$step = step
  }
  if (!missing(k)) {
    args$k = k
  }
  if (!missing(minarea)) {
    args$minarea = minarea
  }
  if (!is.null(centers_arg)) {
    args$centers = centers_arg
  }

  slic_runner = if (iter == 0) sc_slic_points else sc_slic
  slic_sf = do.call(slic_runner, args)
  if (isTRUE(trans$did_transform)) {
    names_x = names(x)
    slic_sf = .supercells_transform_from_lab(slic_sf, names_x)
  }
  return(slic_sf)
}

.supercells_transform_to_lab = function(x, transform) {
  if (is.null(transform)) {
    return(list(x = x, did_transform = FALSE))
  }
  if (!identical(transform, "to_LAB")) {
    stop("The 'transform' argument must be NULL or 'to_LAB'", call. = FALSE)
  }
  if (terra::nlyr(x) < 3) {
    stop("The 'transform = \"to_LAB\"' option requires at least three layers", call. = FALSE)
  }
  if (terra::nlyr(x) > 3) {
    warning("The provided raster has more than three layers: only the first three were kept for calculations", call. = FALSE)
    x = x[[1:3]]
  }
  vals = terra::values(x, mat = TRUE, na.rm = FALSE)
  vals = vals / 255
  vals = grDevices::convertColor(vals, from = "sRGB", to = "Lab")
  x_lab = x
  terra::values(x_lab) = vals
  return(list(x = x_lab, did_transform = TRUE))
}

.supercells_transform_from_lab = function(slic_sf, names_x) {
  value_cols = names_x
  if (is.null(value_cols) || length(value_cols) != 3 || !all(value_cols %in% names(slic_sf))) {
    value_cols = setdiff(names(slic_sf), c("supercells", "x", "y", "geometry"))
    if (length(value_cols) < 3) {
      return(slic_sf)
    }
    value_cols = value_cols[1:3]
  }

  vals = as.matrix(sf::st_drop_geometry(slic_sf)[, value_cols, drop = FALSE])
  storage.mode(vals) = "double"
  vals = grDevices::convertColor(vals, from = "Lab", to = "sRGB") * 255
  vals = pmin(pmax(vals, 0), 255)
  slic_sf[, value_cols] = vals
  return(slic_sf)
}
