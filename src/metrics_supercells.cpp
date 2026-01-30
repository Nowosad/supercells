#include "distances.h"
#include "metrics_helpers.h"
#include "cpp11.hpp"
#include "cpp11/list.hpp"
#include "cpp11/matrix.hpp"
#include <cmath>
#include <vector>

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

  cpp11::writable::doubles mean_value(ncenters);
  cpp11::writable::doubles mean_spatial(ncenters);
  cpp11::writable::doubles mean_combined(ncenters);
  cpp11::writable::doubles balance_ratio(ncenters);
  cpp11::writable::doubles mean_value_scaled(ncenters);

  // Convert sums into per-cluster means and compactness ratios
  for (int i = 0; i < ncenters; i++) {
    if (count[i] > 0) {
      double mv = sum_value[i] / count[i];
      double ms = sum_spatial[i] / count[i];
      double mc = sum_combined[i] / count[i];
      mean_value[i] = mv;
      mean_spatial[i] = ms;
      mean_combined[i] = mc;
      double denom = compactness;
      if (adaptive_compactness) {
        denom = max_value_dist[i];
      }
      mean_value_scaled[i] = (denom != 0.0) ? mv / denom : NA_REAL;
      if (denom != 0.0 && step != 0 && ms > 0.0) {
        balance_ratio[i] = (mv / denom) / (ms / step);
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
