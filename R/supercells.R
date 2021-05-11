#' Creates supercells
#'
#' @param x a
#' @param step a
#' @param nc a
#' @param dist_fun a
#' @param clean a
#' @param iter a
#'
#' @return a
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
