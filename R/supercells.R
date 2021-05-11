#' Creates supercells
#'
#' Creates supercells based on single- or multi-band spatial raster data. It uses a modified version of the SLIC Superpixel algorithm by Achanta et al. (2012), allowing specification of a distance function.
#'
#' @param x An object of class SpatRaster (terra)
#' @param k The number of supercells desired by the user (the output number can be slightly different!)
#' @param compactness A compactness value. Larger values cause clusters to be more compact/even (square).
#' A compactness value depends on the range of input cell values and selected distance measure.
#' @param dist_fun A distance function. Currently implemented distance functions are "euclidean" and "jensen_shannon".
#' @param clean Should connectivity of the supercells be enforced?
#' @param iter The number of iterations performed to create the output.
#'
#' @return An sf object with several columns: (1) supercells - an id of each supercell, (2) y and x coordinates, (3) one or more columns with average values of given variables in each supercell
#'
#' @references Achanta, R., Shaji, A., Smith, K., Lucchi, A., Fua, P., & Süsstrunk, S. (2012). SLIC Superpixels Compared to State-of-the-Art Superpixel Methods. IEEE Transactions on Pattern Analysis and Machine Intelligence, 34(11), 2274–2282. https://doi.org/10.1109/tpami.2012.120
#' @export
#'
#' @examples
#' #a
supercells = function(x, k, compactness, dist_fun = "euclidean", clean = TRUE, iter = 10){
  centers = TRUE
  if (!inherits(x, "SpatRaster")){
    stop("The SpatRaster class is expected as an input")
  }
  mat = dim(x)[1:2]
  mode(mat) = "integer"
  vals = as.matrix(terra::as.data.frame(x, cell = TRUE)[-1])
  slic = run_slic(mat, vals, k = k, nc = compactness, con = clean, centers = centers, type = dist_fun, iter = iter)

  slic_sf = terra::rast(slic[[1]])
  terra::ext(slic_sf) = terra::ext(x)
  terra::crs(slic_sf) = terra::crs(x)
  slic_sf = sf::st_as_sf(terra::as.polygons(slic_sf, dissolve = TRUE))
  # if (centers){
    slic_sf = cbind(slic_sf, slic[[2]])
    names(slic_sf) = c("supercells", "y", "x", "geometry")
    slic_sf[["supercells"]] = slic_sf[["supercells"]] + 1
    slic_sf[["x"]] = as.vector(terra::ext(x))[[1]] + (slic_sf[["x"]] * terra::res(x)[[1]]) + (terra::res(x)[[1]]/2)
    slic_sf[["y"]] = as.vector(terra::ext(x))[[4]] - (slic_sf[["y"]] * terra::res(x)[[2]]) - (terra::res(x)[[1]]/2)
    colnames(slic[[3]]) = names(x)
    slic_sf = cbind(slic_sf, slic[[3]])
  # }
  return(slic_sf)
}

# plot(x)
# plot(slic_sf, add = TRUE, col = NA)
# plot(vect(st_drop_geometry(slic_sf[c("x", "y")]), geom = c("x", "y")), add = TRUE, col = "red")
