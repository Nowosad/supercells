#include "metrics_helpers.h"
#include "cpp11.hpp"
#include "cpp11/list.hpp"
#include "cpp11/matrix.hpp"

[[cpp11::register]]
cpp11::list sc_metrics_global_cpp(cpp11::integers_matrix<> clusters,
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

  double mean_value_sum = 0.0;
  double mean_spatial_sum = 0.0;
  double mean_combined_sum = 0.0;
  double balance_log_sum = 0.0;
  int balance_count = 0;
  double mean_value_scaled_sum = 0.0;
  int active_clusters = 0;

  // Aggregate per-cluster means into global summaries
  for (int i = 0; i < ncenters; i++) {
    if (sums.count[i] <= 0) {
      continue;
    }
    double value_mean = sums.value[i] / sums.count[i];
    double spatial_mean = sums.spatial[i] / sums.count[i];
    double combined_mean = sums.combined[i] / sums.count[i];
    mean_value_sum += value_mean;
    mean_spatial_sum += spatial_mean;
    mean_combined_sum += combined_mean;

    double denominator = data.denominators[i];
    if (denominator != 0.0) {
      mean_value_scaled_sum += value_mean / denominator;
    }
    if (denominator != 0.0 && step != 0 && spatial_mean > 0.0) {
      double ratio = (value_mean / denominator) / (spatial_mean / step);
      if (ratio > 0.0) {
        balance_log_sum += std::log(ratio);
        balance_count += 1;
      }
    }
    active_clusters += 1;
  }

  double mean_value = NA_REAL;
  double mean_spatial = NA_REAL;
  double mean_combined = NA_REAL;
  double balance_log_mean = NA_REAL;
  double mean_value_scaled = NA_REAL;
  if (active_clusters > 0) {
    mean_value = mean_value_sum / active_clusters;
    mean_spatial = mean_spatial_sum / active_clusters;
    mean_combined = mean_combined_sum / active_clusters;
    mean_value_scaled = mean_value_scaled_sum / active_clusters;
  }
  if (balance_count > 0) {
    balance_log_mean = balance_log_sum / balance_count;
  }

  cpp11::writable::list result(6);
  result.names() = {"n_supercells", "mean_value_dist", "mean_spatial_dist",
                    "mean_combined_dist", "balance", "mean_value_dist_scaled"};
  result.at(0) = cpp11::as_sexp(active_clusters);
  result.at(1) = cpp11::as_sexp(mean_value);
  result.at(2) = cpp11::as_sexp(mean_spatial);
  result.at(3) = cpp11::as_sexp(mean_combined);
  result.at(4) = cpp11::as_sexp(balance_log_mean);
  result.at(5) = cpp11::as_sexp(mean_value_scaled);
  return result;
}
