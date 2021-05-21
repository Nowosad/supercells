#' Creates supermotifs
#'
#' Creates supermotifs based on the lsp objects from the motif package. It uses a modified version of the SLIC Superpixel algorithm by Achanta et al. (2012), allowing specification of a distance function.
#'
#' @param x An object of class lsp (output of `motif::lsp_signature()` from the motif package)
#' @param k The number of supercells desired by the user (the output number can be slightly different!)
#' @param compactness A compactness value. Larger values cause clusters to be more compact/even (square).
#' A compactness value depends on the range of input cell values and selected distance measure.
#' @param dist_fun A distance function. Currently implemented distance functions are "euclidean" and "jensen_shannon". Default: "jensen_shannon"
#' @param clean Should connectivity of the supercells be enforced?
#' @param iter The number of iterations performed to create the output.
#'
#' @return An sf object with several columns: (1) supercells - an id of each supercell, (2) y and x coordinates, (3) one or more columns with average values of given variables in each supercell
#'
#' @references Achanta, R., Shaji, A., Smith, K., Lucchi, A., Fua, P., & Süsstrunk, S. (2012). SLIC Superpixels Compared to State-of-the-Art Superpixel Methods. IEEE Transactions on Pattern Analysis and Machine Intelligence, 34(11), 2274–2282. https://doi.org/10.1109/tpami.2012.120
#' @export
#'
#' @examples
#' if (requireNamespace("motif", quietly = TRUE) && requireNamespace("stars", quietly = TRUE)) {
#'   library(motif)
#'   library(stars)
#'   landcover = read_stars(system.file("raster/landcover2015.tif", package = "motif"))
#'   coma_output2 = lsp_signature(landcover, type = "cove", window = 10,
#'                                normalization = "pdf", ordered = FALSE)
#'   slic2 = supermotifs(coma_output2, k = 2000, compactness = 0.1, dist_fun = "jensen_shannon")
#'
#'   plot(landcover, reset = FALSE)
#'   plot(st_geometry(slic2), add = TRUE)
#'
#' }
supermotifs = function(x, k, compactness, dist_fun = "jensen_shannon", clean = TRUE, iter = 10){
  centers = TRUE
  output_sig = motif::lsp_restructure(x)
  output_stars = motif::lsp_add_stars(output_sig)
  mat = dim(output_stars)[2:1]
  vals = as.matrix(as.data.frame(output_stars)[-c(1:4)])

  slic = run_slic(mat, vals, k = k, nc = compactness, con = clean,
                  centers = centers, type = dist_fun, iter = iter)

  slic_sf = terra::rast(slic[[1]])
  terra::NAflag(slic_sf) = -1
  input_ext = sf::st_bbox(output_stars)[c(1, 3, 2, 4)]
  terra::ext(slic_sf) = input_ext
  terra::crs(slic_sf) = sf::st_crs(output_stars)$wkt
  slic_sf = sf::st_as_sf(terra::as.polygons(slic_sf, dissolve = TRUE))
  slic_sf = cbind(slic_sf, stats::na.omit(slic[[2]]))
  names(slic_sf) = c("supercells", "y", "x", "geometry")
  slic_sf[["supercells"]] = slic_sf[["supercells"]] + 1

  res_x = stars::st_dimensions(output_stars)[[1]]$delta
  res_y = stars::st_dimensions(output_stars)[[2]]$delta

  slic_sf[["x"]] = as.vector(input_ext)[[1]] + (slic_sf[["x"]] * res_x) + (res_x/2)
  slic_sf[["y"]] = as.vector(input_ext)[[4]] + (slic_sf[["y"]] * res_y) + (res_y/2)
  # colnames(slic[[3]]) = names(x)
  slic_sf = cbind(slic_sf, stats::na.omit(slic[[3]]))
  return(slic_sf)
}
