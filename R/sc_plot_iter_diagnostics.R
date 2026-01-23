#' Plot iteration diagnostics
#'
#' Plot mean distance across iterations for a supercells run
#'
#' @param x A supercells object with an \code{iter_diagnostics} attribute,
#'   or a diagnostics list containing \code{mean_distance}
#'
#' @return Invisibly returns \code{TRUE} when a plot is created
#'
#' @seealso [`sc_slic()`], [`sc_slic_points()`]
#' @export
#'
#' @examples
#' library(supercells)
#' vol = terra::rast(system.file("raster/volcano.tif", package = "supercells"))
#' vol_sc = sc_slic_points(vol, step = 8, compactness = 1, iter_diagnostics = TRUE)
#' sc_plot_iter_diagnostics(vol_sc)
sc_plot_iter_diagnostics = function(x) {

  iter = attr(x, "iter_diagnostics")

  if (is.null(iter) || is.null(iter$mean_distance) || length(iter$mean_distance) == 0) {
    stop("No iter_diagnostics found", call. = FALSE)
  }

  y = iter$mean_distance
  graphics::plot(seq_along(y), y, type = "b",
                 xlab = "Iteration", ylab = "Mean distance",
                 main = "SLIC iteration diagnostics")
  return(invisible(TRUE))
}
