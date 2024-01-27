#' @export
#'
#' @examples
#' library(supercells)
#' # One variable
#'
#' vol = terra::rast(system.file("raster/volcano.tif", package = "supercells"))
#' run_ce(vol, k = 50)
#'
#' RGB variables
#' ortho = terra::rast(system.file("raster/ortho.tif", package = "supercells"))
#' run_ce(ortho, k = 1000)
# transform is missing
run_ce = function(x, k, dist_fun = "euclidean", step){
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
  # prepare distance function (euclidean is the default)
  if (is.character(dist_fun)){
    if (!(dist_fun %in% c("euclidean", "jsd", "dtw", "dtw2d", philentropy::getDistMethods()))){
      stop("The provided distance function ('dist_fun') does not exist!", call. = FALSE)
    }
    dist_name = dist_fun; dist_fun = function() ""
  } else {
    dist_name = ""
  }
  centers = TRUE
  mat = dim(x)[1:2]; mode(mat) = "integer"
  vals = as.matrix(terra::as.data.frame(x, cells = FALSE, na.rm = FALSE))
  mode(vals) = "double"
  # runs the algorithm
  res = run_compactness_estimation(mat, vals = vals, step = step,
                                   dist_name = dist_name, dist_fun = dist_fun,
                                   input_centers = input_centers)
  return(res)
}
