#' Title
#'
#' @param x
#' @param step
#' @param nc
#' @param dist_fun
#' @param clean
#'
#' @return
#' @export
#'
#' @examples
supercell = function(x, step, nc, dist_fun = "euclidean", clean = TRUE){
  if (!inherits(x, "SpatRaster")){
    stop("The SpatRaster class is expected as an input")
  }
  mat = terra::as.matrix(x[[1]], wide = TRUE)
  mode(mat) = "integer"
  vals = as.matrix(terra::as.data.frame(x, cell = TRUE)[-1])
  slic = run_slic(mat, vals, step = step, nc = nc, con = clean, output_type = TRUE, type = dist_fun)

  slic = terra::rast(slic)
  terra::ext(slic) = terra::ext(x)
  terra::crs(slic) = terra::crs(x)
  slic = sf::st_as_sf(terra::as.polygons(slic, dissolve = TRUE))
  return(slic)
}
