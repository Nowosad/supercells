#' SLIC convergence diagnostics
#'
#' Runs SLIC and returns per-iteration mean combined distance.
#' The output can be plotted directly with [plot()].
#'
#' @inheritParams sc_slic
#'
#' @return A data frame with class `sc_slic_convergence` and columns:
#' \describe{
#'   \item{iter}{Iteration number.}
#'   \item{mean_distance}{Mean combined distance across assigned cells at each iteration.}
#' }
#'
#' @seealso [sc_slic()], [plot()]
#' @export
#'
#' @examples
#' library(supercells)
#' vol = terra::rast(system.file("raster/volcano.tif", package = "supercells"))
#' conv = sc_slic_convergence(vol, step = 8, compactness = 5, iter = 10)
#' plot(conv)
sc_slic_convergence = function(x, step = NULL, compactness, dist_fun = "euclidean",
                               avg_fun = "mean", clean = TRUE, minarea, iter = 10,
                               k = NULL, centers = NULL, verbose = 0) {
  if (iter == 0) {
    stop("iter must be > 0 for convergence diagnostics", call. = FALSE)
  }

  prep = .sc_slic_prep_args(
    x, step, compactness, dist_fun, avg_fun, clean, minarea, iter,
    k, centers, outcomes = "values", chunks = FALSE, verbose = verbose
  )

  res = .sc_run_full_raster(
    x = prep$x,
    step = prep$step,
    compactness = prep$compactness,
    dist_name = prep$funs$dist_name,
    dist_fun = prep$funs$dist_fun,
    adaptive_compactness = prep$adaptive_compactness,
    avg_fun_fun = prep$funs$avg_fun_fun,
    avg_fun_name = prep$funs$avg_fun_name,
    clean = prep$clean,
    iter = prep$iter,
    minarea = prep$minarea,
    input_centers = prep$input_centers,
    iter_diagnostics = TRUE,
    verbose = prep$verbose_cpp
  )

  iter_diag = res$iter_diagnostics
  if (is.null(iter_diag) || is.null(iter_diag$mean_distance) || length(iter_diag$mean_distance) == 0) {
    stop("No convergence diagnostics available", call. = FALSE)
  }

  y = as.numeric(iter_diag$mean_distance)
  out = data.frame(iter = seq_along(y), mean_distance = y)
  class(out) = c("sc_slic_convergence", class(out))
  out
}

#' @export
plot.sc_slic_convergence = function(x, ...) {
  if (!all(c("iter", "mean_distance") %in% names(x))) {
    stop("x must contain 'iter' and 'mean_distance' columns", call. = FALSE)
  }
  graphics::plot(
    x$iter, x$mean_distance, type = "b",
    xlab = "Iteration", ylab = "Mean distance",
    main = "SLIC convergence", ...
  )
  invisible(x)
}
