#' Mark step values as meters
#'
#' Creates a units value in meters for use in `step` arguments.
#' Use plain numerics for cell units, and `use_meters()` for map-distance steps.
#'
#' @param x A single positive numeric value.
#'
#' @return A [units::units] object in meters (`m`).
#'
#' @examples
#' use_meters(100)
#'
#' @export
use_meters = function(x) {
  if (!is.numeric(x) || length(x) != 1 || is.na(x) || x <= 0) {
    stop("The 'x' argument must be a single positive number", call. = FALSE)
  }
  units::set_units(as.numeric(x), "m")
}
