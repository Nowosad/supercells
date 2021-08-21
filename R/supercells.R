#' Creates supercells
#'
#' Creates supercells based on single- or multi-band spatial raster data. It uses a modified version of the SLIC Superpixel algorithm by Achanta et al. (2012), allowing specification of a distance function.
#'
#' @param x An object of class SpatRaster (terra)
#' @param k The number of supercells desired by the user (the output number can be slightly different!).
#' You can use either `k` or `step`.
#' @param compactness A compactness value. Larger values cause clusters to be more compact/even (square).
#' A compactness value depends on the range of input cell values and selected distance measure.
#' @param dist_fun A distance function. Currently implemented distance functions are "euclidean", "jsd", "dtw" (dynamic time warping), name of any distance function from the `philentropy` package (see [philentropy::getDistMethods()]), or any user defined function accepting two vectors and returning one value. Default: "euclidean"
#' @param avg_fun An averaging function - how the values of the supercells' centers are calculated?
#' It accepts any fitting R function (e.g., `base::mean()` or `stats::median()`) or one of internally implemented `"mean"` and `"median"`. Default: `"mean"`
#' @param clean Should connectivity of the supercells be enforced?
#' @param iter The number of iterations performed to create the output.
#' @param minarea Specifies the minimal size of a supercell (in cells). Only works when `clean = TRUE`.
#' By default, when `clean = TRUE`, average area (A) is calculated based on the total number of cells divided by a number of superpixels.
#' Next, the minimal size of a supercell equals to A/(2^2) (A is being right shifted)
#' @param step The distance (number of cells) between initial supercells' centers. You can use either `k` or `step`.
#' @param transform Transformation to be performed on the input. Currently implemented is "to_LAB" allowing to convert RGB raster to a raster in the LAB color space. By default, no transformation is performed.
#'
#' @return An sf object with several columns: (1) supercells - an id of each supercell, (2) y and x coordinates, (3) one or more columns with average values of given variables in each supercell
#'
#' @references Achanta, R., Shaji, A., Smith, K., Lucchi, A., Fua, P., & Süsstrunk, S. (2012). SLIC Superpixels Compared to State-of-the-Art Superpixel Methods. IEEE Transactions on Pattern Analysis and Machine Intelligence, 34(11), 2274–2282. https://doi.org/10.1109/tpami.2012.120
#' @references Nowosad, J. Motif: an open-source R tool for pattern-based spatial analysis. Landscape Ecol (2021). https://doi.org/10.1007/s10980-020-01135-0
#' @export
#'
#' @examples
#' library(supercells)
#' library(terra)
#' library(sf)
#' # One variable
#'
#' vol = rast(system.file("raster/volcano.tif", package = "supercells"))
#' vol_slic1 = supercells(vol, k = 50, compactness = 1)
#' plot(vol)
#' plot(st_geometry(vol_slic1), add = TRUE, lwd = 0.2)
#'
#' # RGB variables
#'
#' ortho = rast(system.file("raster/ortho.tif", package = "supercells"))
#' ortho_slic1 = supercells(ortho, k = 1000, compactness = 10, transform = "to_LAB")
#' plot(ortho)
#' plot(st_geometry(ortho_slic1), add = TRUE)
#'
#' ### RGB variables - colored output
#'
#' rgb_to_hex = function(x){
#'   apply(t(x), 2, function(x) rgb(x[1], x[2], x[3], maxColorValue = 255))
#' }
#' avg_colors = rgb_to_hex(st_drop_geometry(ortho_slic1[4:6]))
#'
#' plot(ortho)
#' plot(st_geometry(ortho_slic1), add = TRUE, col = avg_colors)
supercells = function(x, k, compactness, dist_fun = "euclidean", avg_fun = "mean", clean = TRUE,
                      iter = 10, transform = NULL, step, minarea){
  centers = TRUE
  if (!inherits(x, "SpatRaster")){
    stop("The SpatRaster class is expected as an input", call. = FALSE)
  }
  mat = dim(x)[1:2]
  mode(mat) = "integer"
  vals = as.matrix(terra::as.data.frame(x, cell = FALSE, na.rm = FALSE))
  if (!missing(step) && !missing(k)){
    stop("You can specify either k or step, not both", call. = FALSE)
  } else if (missing(step) && missing(k)){
    stop("You need to specify either k or step", call. = FALSE)
  } else if (missing(step)){
    superpixelsize = round((mat[1] * mat[2]) / k + 0.5)
    step = round(sqrt(superpixelsize) + 0.5)
  }
  if (!missing(transform)){
    if (transform == "to_LAB"){
      vals = vals / 255
      vals = grDevices::convertColor(vals, from = "sRGB", to = "Lab")
    }
  }
  if (is.character(avg_fun)){
    avg_fun_name = avg_fun
    avg_fun_fun = function() ""
  } else {
    avg_fun_name = ""
    avg_fun_fun = avg_fun
  }
  if (is.character(dist_fun)){
    if (!(dist_fun %in% c("euclidean", "jsd", "dtw", philentropy::getDistMethods()))){
      stop("The provided distance function ('dist_fun') does not exist!", call. = FALSE)
    }
    dist_type = dist_fun
    dist_fun = function() ""
  } else {
    dist_type = ""
  }
  if (missing(minarea)){
    minarea = 0
  }
  slic = run_slic(mat, vals = vals, step = step, nc = compactness, con = clean,
                  centers = centers, type = dist_type, type_fun = dist_fun,
                  avg_fun_fun = avg_fun_fun, avg_fun_name = avg_fun_name,
                  iter = iter, lims = minarea)
  if (!missing(transform)){
    if (transform == "to_LAB"){
      slic[[3]] = grDevices::convertColor(slic[[3]], from = "Lab", to = "sRGB") * 255
    }
  }
  slic_sf = terra::rast(slic[[1]])
  terra::NAflag(slic_sf) = -1
  terra::ext(slic_sf) = terra::ext(x)
  terra::crs(slic_sf) = terra::crs(x)
  slic_sf = sf::st_as_sf(terra::as.polygons(slic_sf, dissolve = TRUE))
  # if (centers){
    slic_sf = cbind(slic_sf, stats::na.omit(slic[[2]]))
    names(slic_sf) = c("supercells", "x", "y", "geometry")
    slic_sf[["supercells"]] = slic_sf[["supercells"]] + 1
    slic_sf[["x"]] = as.vector(terra::ext(x))[[1]] + (slic_sf[["x"]] * terra::res(x)[[1]]) + (terra::res(x)[[1]]/2)
    slic_sf[["y"]] = as.vector(terra::ext(x))[[4]] - (slic_sf[["y"]] * terra::res(x)[[2]]) - (terra::res(x)[[1]]/2)
    colnames(slic[[3]]) = names(x)
    slic_sf = cbind(slic_sf, stats::na.omit(slic[[3]]))
    slic_sf = suppressWarnings(sf::st_collection_extract(slic_sf, "POLYGON"))

  # }
  return(slic_sf)
}

# plot(x)
# plot(slic_sf, add = TRUE, col = NA)
# plot(vect(st_drop_geometry(slic_sf[c("x", "y")]), geom = c("x", "y")), add = TRUE, col = "red")
# rgb_to_lab = function(x){
#   new_vals = values(logo) / 255
#   new_vals = convertColor(new_vals, from = "sRGB", to = "Lab")
#   values(x) = new_vals
#   x
# }
