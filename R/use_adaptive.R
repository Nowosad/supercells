#' Use adaptive compactness mode
#'
#' Creates a compactness mode object for adaptive compactness.
#' The `"local_max"` method scales value distances separately for each center
#' using the largest feature-space distance within its local search window.
#' The scale is updated before each assignment iteration, while spatial
#' distances remain scaled by `step`.
#'
#' @param method Adaptive compactness method. Currently only `"local_max"` is
#' supported.
#'
#' @return An adaptive compactness mode object for `compactness` arguments.
#'
#' @examples
#' use_adaptive()
#'
#' @export
use_adaptive = function(method = "local_max") {
  if (!is.character(method) || length(method) != 1 || is.na(method) || method != "local_max") {
    stop("The 'method' argument must be 'local_max'", call. = FALSE)
  }
  structure(list(method = method), class = "sc_adaptive")
}
