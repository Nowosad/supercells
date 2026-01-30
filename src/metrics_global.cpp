#include "distances.h"
#include "metrics_helpers.h"
#include "cpp11.hpp"
#include "cpp11/list.hpp"
#include "cpp11/matrix.hpp"
#include <cmath>
#include <vector>

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
  int bands = vals.ncol();
  int ncenters = centers_vals.nrow();

  // Copy center values into std::vector for distance computation
  std::vector<std::vector<double>> centers_vals_vec(ncenters, std::vector<double>(bands));
  for (int i = 0; i < ncenters; i++) {
    for (int nval = 0; nval < bands; nval++) {
      centers_vals_vec[i][nval] = centers_vals(i, nval);
    }
  }

  std::vector<double> sum_value(ncenters, 0.0);
  std::vector<double> sum_spatial(ncenters, 0.0);
  std::vector<double> sum_combined(ncenters, 0.0);
  std::vector<int> count(ncenters, 0);

  std::vector<double> pixel_values;
  pixel_values.reserve(bands);

  std::vector<double> max_value_dist;
  if (adaptive_compactness) {
    max_value_dist = sc_compute_max_value_dist(centers_vals_vec, centers_xy, vals,
                                               rows, cols, bands, step, dist_name, dist_fun);
  }

  // Per-pixel pass: accumulate distances into per-cluster sums
  for (int i = 0; i < cols; i++) {
    for (int j = 0; j < rows; j++) {
      int cid = clusters(j, i);
      if (cid < 0 || cid >= ncenters) {
        continue;
      }

      int ncell = i + (j * cols);
      pixel_values.clear();
      bool has_na = false;
      for (int nval = 0; nval < bands; nval++) {
        double val = vals(ncell, nval);
        pixel_values.push_back(val);
        if (std::isnan(val)) {
          has_na = true;
        }
      }
      if (has_na) {
        continue;
      }

      // Compute value- and spatial-space distances to the cluster center
      double value_dist = get_vals_dist(centers_vals_vec[cid], pixel_values, dist_name, dist_fun);
      double center_x = centers_xy(cid, 0);
      double center_y = centers_xy(cid, 1);
      double y_dist = center_y - j;
      double x_dist = center_x - i;
      double spatial_dist = sqrt((y_dist * y_dist) + (x_dist * x_dist));

      double denom = compactness;
      if (adaptive_compactness) {
        denom = max_value_dist[cid];
      }

      double combined_dist = NA_REAL;
      if (denom != 0.0 && step != 0) {
        double dist1 = value_dist / denom;
        double dist2 = spatial_dist / step;
        combined_dist = sqrt((dist1 * dist1) + (dist2 * dist2));
      }

      // Accumulate sums; combined is only counted when well-defined
      sum_value[cid] += value_dist;
      sum_spatial[cid] += spatial_dist;
      if (!std::isnan(combined_dist)) {
        sum_combined[cid] += combined_dist;
      }
      count[cid] += 1;
    }
  }

  double mean_value_sum = 0.0;
  double mean_spatial_sum = 0.0;
  double mean_combined_sum = 0.0;
  double balance_ratio_sum = 0.0;
  double mean_value_scaled_sum = 0.0;
  int active_clusters = 0;

  // Aggregate per-cluster means into global summaries
  for (int i = 0; i < ncenters; i++) {
    if (count[i] <= 0) {
      continue;
    }
    double mv = sum_value[i] / count[i];
    double ms = sum_spatial[i] / count[i];
    double mc = sum_combined[i] / count[i];
    mean_value_sum += mv;
    mean_spatial_sum += ms;
    mean_combined_sum += mc;

    double denom = compactness;
    if (adaptive_compactness) {
      denom = max_value_dist[i];
    }
    if (denom != 0.0) {
      mean_value_scaled_sum += mv / denom;
    }
    if (denom != 0.0 && step != 0 && ms > 0.0) {
      balance_ratio_sum += (mv / denom) / (ms / step);
    }
    active_clusters += 1;
  }

  double mean_value = NA_REAL;
  double mean_spatial = NA_REAL;
  double mean_combined = NA_REAL;
  double balance_ratio_mean = NA_REAL;
  double mean_value_scaled = NA_REAL;
  if (active_clusters > 0) {
    mean_value = mean_value_sum / active_clusters;
    mean_spatial = mean_spatial_sum / active_clusters;
    mean_combined = mean_combined_sum / active_clusters;
    balance_ratio_mean = balance_ratio_sum / active_clusters;
    mean_value_scaled = mean_value_scaled_sum / active_clusters;
  }

  cpp11::writable::list result(6);
  result.names() = {"n_supercells", "mean_value_dist", "mean_spatial_dist",
                    "mean_combined_dist", "balance", "mean_value_dist_scaled"};
  result.at(0) = cpp11::as_sexp(active_clusters);
  result.at(1) = cpp11::as_sexp(mean_value);
  result.at(2) = cpp11::as_sexp(mean_spatial);
  result.at(3) = cpp11::as_sexp(mean_combined);
  result.at(4) = cpp11::as_sexp(balance_ratio_mean);
  result.at(5) = cpp11::as_sexp(mean_value_scaled);
  return result;
}
