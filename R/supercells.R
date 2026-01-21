#' Creates supercells
#'
#' Creates supercells based on single- or multi-band spatial raster data. It uses a modified version of the SLIC Superpixel algorithm by Achanta et al. (2012), allowing specification of a distance function.
#'
#' @param x An object of class SpatRaster (terra) or class stars (stars)
#' @param k The number of supercells desired by the user (the output number can be slightly different!).
#' You can use either `k` or `step`.
#' It is also possible to provide a set of points (an `sf` object) as `k` together with the `step` value to create custom cluster centers.
#' @param compactness A compactness value. Larger values cause clusters to be more compact/even (square).
#' A compactness value depends on the range of input cell values and selected distance measure.
#' @param dist_fun A distance function. Currently implemented distance functions are `"euclidean"`, `"jsd"`, `"dtw"` (dynamic time warping), name of any distance function from the `philentropy` package (see [philentropy::getDistMethods()]; "log2" is used in this case), or any user defined function accepting two vectors and returning one value. Default: `"euclidean"`
#' @param avg_fun An averaging function - how the values of the supercells' centers are calculated? The algorithm internally implements common functions `"mean"` and `"median"` (provided with quotation marks), but also accepts any fitting R function (e.g., `base::mean()` or `stats::median()`, provided as plain function name: `mean`). Default: `"mean"`. See details for more information.
#' @param clean Should connectivity of the supercells be enforced?
#' @param iter The number of iterations performed to create the output.
#' @param minarea Specifies the minimal size of a supercell (in cells). Only works when `clean = TRUE`.
#' By default, when `clean = TRUE`, average area (A) is calculated based on the total number of cells divided by a number of supercells
#' Next, the minimal size of a supercell equals to A/(2^2) (A is being right shifted)
#' @param step The distance (number of cells) between initial supercells' centers. You can use either `k` or `step`.
#' @param transform Transformation to be performed on the input. By default, no transformation is performed. Currently available transformation is "to_LAB": first, the conversion from RGB to the LAB color space is applied, then the supercells algorithm is run, and afterward, a reverse transformation is performed on the obtained results. (This argument is experimental and may be removed in the future).
#' @param metadata Logical. If `TRUE`, the output object will have metadata columns ("supercells", "x", "y"). If `FALSE`, the output object will not have metadata columns.
#' @param chunks Should the input (`x`) be split into chunks before deriving supercells? Either `FALSE` (default), `TRUE` (only large input objects are split), or a numeric value (representing the side length of the chunk in the number of cells).
#' @param future Should the future package be used for parallelization of the calculations? Default: `FALSE`. If `TRUE`, you also need to specify `future::plan()`.
#' @param verbose An integer specifying the level of text messages printed during calculations. 0 means no messages (default), 1 provides basic messages (e.g., calculation stage).
#'
#' @details
#' If you want to use additional arguments for the averaging function (`avg_fun`), you can create a custom function. For example, if you want to calculate the mean by removing missing values, you can use the following code: `my_mean = function(x) mean(x, na.rm = TRUE)` and then provide `avg_fun = my_mean.`
#'
#' @return An sf object with several columns: (1) supercells - an id of each supercell, (2) y and x coordinates, (3) one or more columns with average values of given variables in each supercell.
#'
#' @references Achanta, R., Shaji, A., Smith, K., Lucchi, A., Fua, P., & Süsstrunk, S. (2012). SLIC Superpixels Compared to State-of-the-Art Superpixel Methods. IEEE Transactions on Pattern Analysis and Machine Intelligence, 34(11), 2274–2282. https://doi.org/10.1109/tpami.2012.120
#' @references Nowosad, J. Motif: an open-source R tool for pattern-based spatial analysis. Landscape Ecol (2021). https://doi.org/10.1007/s10980-020-01135-0
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
                      iter = 10, transform = NULL, step, minarea, metadata = TRUE, chunks = FALSE, future = FALSE, verbose = 0){
  if (!inherits(x, "SpatRaster")){
    if (inherits(x, "stars")){
      x = terra::rast(x)
    } else{
      stop("The SpatRaster class is expected as an input", call. = FALSE)
    }
  }
  # prepare initial supercells' centers
  input_centers = matrix(c(0L, 0L), ncol = 2)
  if (!missing(k) && inherits(k, "sf")){
    if (chunks > 0){
      stop(call. = FALSE, "Chunks cannot be used for custom cluster centers!")
    }
    input_centers = centers_to_dims(x, k)
  } else if (!missing(step) && !missing(k)){
    stop("You can specify either k or step, not both", call. = FALSE)
  } else if (missing(step) && missing(k)){
    stop("You need to specify either k or step", call. = FALSE)
  } else if (missing(step)){
    mat = dim(x)[1:2]; mode(mat) = "integer"
    superpixelsize = round((mat[1] * mat[2]) / k + 0.5)
    step = round(sqrt(superpixelsize) + 0.5)
  }
  # prepare averaging function (mean is the default)
  if (is.character(avg_fun)){
    avg_fun_name = avg_fun; avg_fun_fun = function() ""
  } else {
    avg_fun_name = "";      avg_fun_fun = avg_fun
  }
  # prepare distance function (euclidean is the default)
  if (is.character(dist_fun)){
    if (!(dist_fun %in% c("euclidean", "jsd", "dtw", "dtw2d", philentropy::getDistMethods()))){
      stop("The provided distance function ('dist_fun') does not exist!", call. = FALSE)
    }
    dist_name = dist_fun; dist_fun = function() ""
  } else {
    dist_name = ""
  }
  # prepare minarea
  if (missing(minarea)){
    minarea = 0
  } else if (minarea > step^2) {
    warning("The provided minarea value is larger than than the average supercell (step^2). The connectivity cleaning is likely to fail.", call. = FALSE)
  }
  # disables cleaning if iter = 0
  if (iter == 0){
    clean = FALSE
  }
  # get extents of chunks
  chunk_ext = prep_chunks_ext(dim(x), limit = chunks)
  # run the algorithm on chunks
  if (future){
    if (in_memory(x)){
        names_x = names(x)
        x = terra::writeRaster(x, tempfile(fileext = ".tif"))
        names(x) = names_x
    }
    if (!in_memory(x)){
      x = terra::sources(x)[[1]]
    }
    oopts = options(future.globals.maxSize = +Inf)
    on.exit(options(oopts))
    slic_sf = future.apply::future_apply(chunk_ext, MARGIN = 1, run_slic_chunks, x = x,
                                         step = step, compactness = compactness, dist_name = dist_name,
                                         dist_fun = dist_fun, avg_fun_fun = avg_fun_fun, avg_fun_name = avg_fun_name,
                                         clean = clean, iter = iter, minarea = minarea, transform = transform,
                                         input_centers = input_centers, verbose = verbose,
                                         future.seed = TRUE)
  } else {
    slic_sf = apply(chunk_ext, MARGIN = 1, run_slic_chunks, x = x,
                                         step = step, compactness = compactness, dist_name = dist_name,
                                         dist_fun = dist_fun, avg_fun_fun = avg_fun_fun, avg_fun_name = avg_fun_name,
                                         clean = clean, iter = iter, minarea = minarea, transform = transform,
                                         input_centers = input_centers, verbose = verbose)
  }
  # combines the chunks results by updating supercells ids
  slic_sf = update_supercells_ids(slic_sf)
  # removes metadata columns if metadata = FALSE
  if (isFALSE(metadata)){
    slic_sf = slic_sf[, -which(names(slic_sf) %in% c("supercells", "x", "y"))]
  }
  # returns the result
  return(slic_sf)
}
