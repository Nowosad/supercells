#' Title
#'
#' @param x a
#' @param step a
#' @param nc a
#' @param dist_fun a
#' @param clean a
#'
#' @return a
#' @export
#'
#' @examples
#' #a
supercell = function(x, step, nc, dist_fun = "euclidean", clean = TRUE, centers = FALSE, iter = 10){
  if (!inherits(x, "SpatRaster")){
    stop("The SpatRaster class is expected as an input")
  }
  mat = dim(x)[1:2]
  mode(mat) = "integer"
  vals = as.matrix(terra::as.data.frame(x, cell = TRUE)[-1])
  slic = run_slic(mat, vals, step = step, nc = nc, con = clean, centers = centers, type = dist_fun, iter = iter)

  slic_sf = terra::rast(slic[[1]])
  terra::ext(slic_sf) = terra::ext(x)
  terra::crs(slic_sf) = terra::crs(x)
  slic_sf = sf::st_as_sf(terra::as.polygons(slic_sf, dissolve = TRUE))
  slic_sf = cbind(slic_sf, slic[[2]])
  names(slic_sf) = c("supercell", "x", "y", "geometry")
  slic_sf[["x"]] = as.vector(terra::ext(x))[[1]] + (slic_sf[["x"]] * terra::res(x)[[1]])
  slic_sf[["y"]] = as.vector(terra::ext(x))[[3]] + (slic_sf[["y"]] * terra::res(x)[[2]])
  return(slic_sf)
}



#' @export
supermotif = function(x, step, nc, dist_fun = "euclidean", clean = TRUE){
  output_sig = lsp_restructure(x)
  output_stars = motif::lsp_add_stars(output_sig)
  mat = dim(output_stars)[2:1]
  vals = as.matrix(as.data.frame(output_stars)[-c(1:4)])

  slic = run_slic(mat, vals, step = step, nc = nc, con = clean, output_type = TRUE, type = dist_fun)

  slic = terra::rast(slic)
  terra::ext(slic) = sf::st_bbox(output_stars)[c(1, 3, 2, 4)]
  terra::crs(slic) = sf::st_crs(output_stars)$wkt
  slic = sf::st_as_sf(terra::as.polygons(slic, dissolve = TRUE))
  return(slic)
}

lsp_restructure = function(x){
  x_attr = attributes(x)
  nc = ncol(x$signature[[1]])

  unnested_signature = matrix(unlist(x$signature, use.names = FALSE),
                              ncol = nc, byrow = TRUE)
  colnames(unnested_signature) = paste0("X", seq_len(nc))
  unnested_signature = tibble::as_tibble(unnested_signature)

  x["signature"] = NULL

  x = tibble::as_tibble(cbind(x, unnested_signature))
  x_attr$names = names(x)
  attributes(x) = x_attr

  x
}
# a = supercell(rast(landscape), 3, 1)
# plot(rast(landscape))
# plot(a, add = TRUE, col = NA)
# format(object.size(as.matrix(dog2[[1]], wide = T)), "MB")
# format(object.size(as.matrix(terra::as.data.frame(dog2, cell = TRUE)[-1])), "MB")

