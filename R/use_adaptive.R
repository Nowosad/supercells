#' Use adaptive compactness mode
#'
#' Creates a compactness mode object for adaptive compactness.
#' The `"local_max"` method corresponds to SLIC0-style local scaling,
#' where compactness is adapted using local maximum value distances.
#'
#' @param method Adaptive compactness method. Currently only `"local_max"` is supported
#' (SLIC0-style).
#'
#' @return An adaptive compactness mode object for `compactness` arguments.
#'
#' @examples
#' use_adaptive()
#'
#' @export
use_adaptive = function(method = "local_max") {
  if (!is.character(method) || length(method) != 1 || is.na(method) || method != "local_max") {
    stop("The 'method' argument must be 'local_max' (SLIC0-style)", call. = FALSE)
  }
  structure(list(method = method), class = "sc_adaptive")
}
