#include "distances.h"
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
                                  std::string dist_name, cpp11::function dist_fun) {
  int rows = clusters.nrow();
  int cols = clusters.ncol();
  int bands = vals.ncol();
  int ncenters = centers_vals.nrow();

  std::vector<std::vector<double>> centers_vals_vec(ncenters, std::vector<double>(bands));
  for (int i = 0; i < ncenters; i++) {
    for (int nval = 0; nval < bands; nval++) {
      centers_vals_vec[i][nval] = centers_vals(i, nval);
    }
  }

  auto dist_cb = [&](const std::vector<double>& values1,
                     const std::vector<double>& values2) -> double {
    return get_vals_dist(values1, values2, dist_name, dist_fun);
  };

  std::vector<double> sum_value(ncenters, 0.0);
  std::vector<double> sum_spatial(ncenters, 0.0);
  std::vector<double> sum_combined(ncenters, 0.0);
  std::vector<int> count(ncenters, 0);

  std::vector<double> pixel_values;
  pixel_values.reserve(bands);

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

      double value_dist = dist_cb(centers_vals_vec[cid], pixel_values);
      double center_x = centers_xy(cid, 0);
      double center_y = centers_xy(cid, 1);
      double y_dist = center_y - j;
      double x_dist = center_x - i;
      double spatial_dist = sqrt((y_dist * y_dist) + (x_dist * x_dist));

      double combined_dist = NA_REAL;
      if (compactness != 0.0 && step != 0) {
        double dist1 = value_dist / compactness;
        double dist2 = spatial_dist / step;
        combined_dist = sqrt((dist1 * dist1) + (dist2 * dist2));
      }

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
  double compact_ratio_sum = 0.0;
  int active_clusters = 0;

  double w_value_sum = 0.0;
  double w_spatial_sum = 0.0;
  double w_combined_sum = 0.0;
  int total_count = 0;

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

    if (compactness != 0.0 && step != 0 && ms > 0.0) {
      compact_ratio_sum += (mv / compactness) / (ms / step);
    }

    w_value_sum += mv * count[i];
    w_spatial_sum += ms * count[i];
    w_combined_sum += mc * count[i];
    total_count += count[i];
    active_clusters += 1;
  }

  double mean_value = NA_REAL;
  double mean_spatial = NA_REAL;
  double mean_combined = NA_REAL;
  double compact_ratio_mean = NA_REAL;
  double mean_value_w = NA_REAL;
  double mean_spatial_w = NA_REAL;
  double mean_combined_w = NA_REAL;

  if (active_clusters > 0) {
    mean_value = mean_value_sum / active_clusters;
    mean_spatial = mean_spatial_sum / active_clusters;
    mean_combined = mean_combined_sum / active_clusters;
    compact_ratio_mean = compact_ratio_sum / active_clusters;
  }
  if (total_count > 0) {
    mean_value_w = w_value_sum / total_count;
    mean_spatial_w = w_spatial_sum / total_count;
    mean_combined_w = w_combined_sum / total_count;
  }

  cpp11::writable::list result(8);
  result.names() = {"n_supercells", "mean_value_dist", "mean_spatial_dist",
                    "mean_combined_dist", "compactness_ratio_mean",
                    "mean_value_dist_w", "mean_spatial_dist_w", "mean_combined_dist_w"};
  result.at(0) = cpp11::as_sexp(active_clusters);
  result.at(1) = cpp11::as_sexp(mean_value);
  result.at(2) = cpp11::as_sexp(mean_spatial);
  result.at(3) = cpp11::as_sexp(mean_combined);
  result.at(4) = cpp11::as_sexp(compact_ratio_mean);
  result.at(5) = cpp11::as_sexp(mean_value_w);
  result.at(6) = cpp11::as_sexp(mean_spatial_w);
  result.at(7) = cpp11::as_sexp(mean_combined_w);
  return result;
}
