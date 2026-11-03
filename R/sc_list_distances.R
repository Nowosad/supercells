# Canonical metadata for supported named distance measures.
.sc_distance_registry = utils::read.delim(text = '
name;label;input;compares;effect_of_rescaling;guidance;tuning_scale
euclidean;Euclidean;numeric values;root-sum-square layer differences;changes distance;Scale layers to comparable units or ranges.;sqrt_layers
squared_euclidean;Squared Euclidean;numeric values;summed squared layer differences;changes distance;Scale layers to comparable units or ranges.;layers
manhattan;Manhattan;numeric values;summed absolute layer differences;changes distance;Scale layers to comparable units or ranges.;layers
chebyshev;Chebyshev;numeric values;largest absolute layer difference;changes distance;Scale layers to comparable units or ranges.;none
canberra;Canberra;numeric values;relative layer differences;changes distance;Layer scaling has no effect.;layers
lorentzian;Lorentzian;numeric values;summed log-scaled layer differences;changes distance;Scale layers to comparable units or ranges.;layers

soergel;Soergel;nonnegative values;difference relative to larger values;changes distance;Scale layers to set their influence.;none
wavehedges;Wave Hedges;nonnegative values;summed relative layer differences;changes distance;Layer scaling has no effect.;layers
czekanowski;Czekanowski;nonnegative values;difference relative to total magnitude;changes distance;Scale layers to set their influence.;none

squared_chord;Squared chord;nonnegative values;differences after square-root transformation;changes distance;Scale layers to set their influence.;layers
squared_chi;Squared chi-square;nonnegative values;squared differences relative to combined values;changes distance;Scale layers to set their influence.;layers
clark;Clark;nonnegative values;combined relative layer differences;changes distance;Layer scaling has no effect.;sqrt_layers
jsd;Jensen-Shannon divergence;probabilities summing to 1;difference between probability profiles;invalid unless values sum to 1;Normalize every cell to sum to 1 and use mean centres.;none
non-intersection;Non-intersection;probabilities summing to 1;probability assigned differently across layers;invalid unless values sum to 1;Normalize every cell to sum to 1 and use mean centres.;none

cosine_distance;Cosine distance;numeric values;pattern across layers, ignoring overall scale;no effect;Scale layers only to change their influence.;none

dtw;Dynamic time warping;ordered numeric sequence;sequence difference after alignment;changes distance;Arrange layers in sequence order.;layers
dtw2d;Two-dimensional dynamic time warping;two equal-length coordinate sequences;2D trajectory difference after alignment;changes distance;Arrange layers by dimension, with sequence order within each dimension.;layers
', sep = ";", quote = "", strip.white = TRUE)

#' List supported distance measures
#'
#' Returns information about the distance measures that can be provided by name
#' to the `dist_fun` argument. The table describes their input
#' requirements, what they compare, and practical guidance.
#'
#' Other measures can be supplied explicitly as a function to `dist_fun`.
#' In that case, the function must take two numeric vectors and return one
#' finite, nonnegative value, with smaller values meaning greater similarity.
#'
#' For `"cosine_distance"`, two zero profiles have distance zero. A zero
#' profile and a nonzero profile have distance one.
#'
#' @return A data frame with one row per supported measure and columns:
#' \describe{
#'   \item{`name`}{Value accepted by `dist_fun`.}
#'   \item{`label`}{Measure name.}
#'   \item{`input`}{Type of the input vector.}
#'   \item{`compares`}{What aspect of the vectors the measure compares.}
#'   \item{`effect_of_rescaling`}{Effect of rescaling the input values.}
#'   \item{`guidance`}{Practical usage guidance.}
#' }
#'
#' @examples
#' sc_list_distances()
#'
#' @importFrom philentropy dist_one_one
#' @export
sc_list_distances = function() {
  .sc_distance_registry[names(.sc_distance_registry) != "tuning_scale"]
}
