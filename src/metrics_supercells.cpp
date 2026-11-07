#include "metrics_helpers.h"
#include "cpp11.hpp"
#include "cpp11/list.hpp"
#include "cpp11/matrix.hpp"

[[cpp11::register]]
cpp11::list sc_metrics_supercells_cpp(cpp11::integers_matrix<> clusters,
                                      cpp11::doubles_matrix<> centers_xy,
                                      cpp11::doubles_matrix<> centers_vals,
                                      cpp11::doubles_matrix<> vals,
                                      int step, double compactness,
                                      bool adaptive_compactness,
                                      std::string dist_name, cpp11::function dist_fun) {
  int rows = clusters.nrow();
  int cols = clusters.ncol();
  int ncenters = centers_vals.nrow();

  ScMetricData data = sc_prepare_metric_data(
    centers_xy, centers_vals, vals, rows, cols, step, compactness,
    adaptive_compactness, dist_name, dist_fun
  );
  ScMetricSums sums = sc_accumulate_metrics(
    clusters, centers_xy, vals, data, step, dist_name, dist_fun
  );

  cpp11::writable::doubles mean_value(ncenters);
  cpp11::writable::doubles mean_spatial(ncenters);
  cpp11::writable::doubles mean_combined(ncenters);
  cpp11::writable::doubles balance_ratio(ncenters);
  cpp11::writable::doubles mean_value_scaled(ncenters);

  // Convert sums into per-cluster means and compactness ratios
  for (int i = 0; i < ncenters; i++) {
    if (sums.count[i] > 0) {
      double value_mean = sums.value[i] / sums.count[i];
      double spatial_mean = sums.spatial[i] / sums.count[i];
      double combined_mean = sums.combined[i] / sums.count[i];
      mean_value[i] = value_mean;
      mean_spatial[i] = spatial_mean;
      mean_combined[i] = combined_mean;
      double denominator = data.denominators[i];
      mean_value_scaled[i] = (denominator != 0.0) ? value_mean / denominator : NA_REAL;
      if (denominator != 0.0 && step != 0 && spatial_mean > 0.0) {
        balance_ratio[i] = (value_mean / denominator) / (spatial_mean / step);
      } else {
        balance_ratio[i] = NA_REAL;
      }

    } else {
      mean_value[i] = NA_REAL;
      mean_spatial[i] = NA_REAL;
      mean_combined[i] = NA_REAL;
      balance_ratio[i] = NA_REAL;
      mean_value_scaled[i] = NA_REAL;
    }
  }

  cpp11::writable::list result(5);
  result.names() = {"mean_value_dist", "mean_spatial_dist", "mean_combined_dist",
                    "balance", "mean_value_dist_scaled"};
  result.at(0) = mean_value;
  result.at(1) = mean_spatial;
  result.at(2) = mean_combined;
  result.at(3) = balance_ratio;
  result.at(4) = mean_value_scaled;
  return result;
}
