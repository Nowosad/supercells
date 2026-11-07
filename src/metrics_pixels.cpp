#include "metrics_helpers.h"
#include "cpp11.hpp"
#include "cpp11/list.hpp"
#include "cpp11/matrix.hpp"

[[cpp11::register]]
cpp11::list sc_metrics_pixels_cpp(cpp11::integers_matrix<> clusters,
                                  cpp11::doubles_matrix<> centers_xy,
                                  cpp11::doubles_matrix<> centers_vals,
                                  cpp11::doubles_matrix<> vals,
                                  int step, double compactness,
                                  bool adaptive_compactness,
                                  std::string dist_name, cpp11::function dist_fun) {
  int rows = clusters.nrow();
  int cols = clusters.ncol();

  cpp11::writable::doubles_matrix<> per_pixel_spatial(rows, cols);
  cpp11::writable::doubles_matrix<> per_pixel_value(rows, cols);
  cpp11::writable::doubles_matrix<> per_pixel_combined(rows, cols);
  cpp11::writable::doubles_matrix<> per_pixel_value_scaled(rows, cols);

  ScMetricData data = sc_prepare_metric_data(
    centers_xy, centers_vals, vals, rows, cols, step, compactness,
    adaptive_compactness, dist_name, dist_fun
  );
  std::vector<double> pixel_values;
  pixel_values.reserve(vals.ncol());

  for (int col = 0; col < cols; col++) {
    for (int row = 0; row < rows; row++) {
      int cid = clusters(row, col);
      ScPixelMetrics metrics = {NA_REAL, NA_REAL, NA_REAL, NA_REAL};
      sc_compute_pixel_metrics(
        cid, col, row, cols, step, centers_xy, vals, data,
        dist_name, dist_fun, pixel_values, metrics
      );
      per_pixel_spatial(row, col) = metrics.spatial;
      per_pixel_value(row, col) = metrics.value;
      per_pixel_combined(row, col) = metrics.combined;
      per_pixel_value_scaled(row, col) = metrics.value_scaled;
    }
  }

  cpp11::writable::list result(4);
  result.names() = {"spatial", "value", "combined", "value_scaled"};
  result.at(0) = per_pixel_spatial;
  result.at(1) = per_pixel_value;
  result.at(2) = per_pixel_combined;
  result.at(3) = per_pixel_value_scaled;
  return result;
}
