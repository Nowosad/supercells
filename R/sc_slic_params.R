#' Get stored `sc_slic()` parameters
#'
#' Returns key `sc_slic()` parameters stored as attributes on a supercells object.
#'
#' @param sc An sf object returned by [sc_slic()].
#'
#' @return A one-row data.frame with columns:
#' `step`, `compactness`, `compactness_method`, and `dist_fun`.
#' The `dist_fun` column is character; custom distance functions are stored as
#' `NA`.
#'
#' @seealso [sc_slic()], [sc_slic_set_params()]
#' @export
sc_slic_get_params = function(sc) {
  dist_fun = attr(sc, "dist_fun")
  if (is.function(dist_fun)) {
    dist_fun = NA_character_
  }
  data.frame(
    step = attr(sc, "step"),
    compactness = attr(sc, "compactness"),
    compactness_method = attr(sc, "compactness_method"),
    dist_fun = dist_fun
  )
}

#' Set stored `sc_slic()` parameters
#'
#' Writes key `sc_slic()` parameters to attributes on a supercells object.
#'
#' @param sc An sf object.
#' @param params A data.frame, typically from [sc_slic_get_params()]. Only the
#' first row is used.
#'
#' @return The input object with updated attributes.
#'
#' @seealso [sc_slic()], [sc_slic_get_params()]
#' @export
sc_slic_set_params = function(sc, params) {
  expected_cols = c("step", "compactness", "compactness_method", "dist_fun")
  if (!all(expected_cols %in% names(params))) {
    stop("params must contain columns: step, compactness, compactness_method, dist_fun", call. = FALSE)
  }

  attr(sc, "step") = params[["step"]][1]
  attr(sc, "compactness") = params[["compactness"]][1]
  attr(sc, "compactness_method") = params[["compactness_method"]][1]
  dist_fun = params[["dist_fun"]][1]
  if (is.na(dist_fun)) {
    attr(sc, "dist_fun") = NULL
  } else {
    attr(sc, "dist_fun") = dist_fun
  }

  sc
}
