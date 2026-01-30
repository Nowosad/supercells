#include "distances.h"
#include "metrics_helpers.h"
#include "cpp11.hpp"
#include "cpp11/list.hpp"
#include "cpp11/matrix.hpp"
#include <cmath>
#include <vector>

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
  int bands = vals.ncol();
  int ncenters = centers_vals.nrow();

  cpp11::writable::doubles_matrix<> per_pixel_spatial(rows, cols);
  cpp11::writable::doubles_matrix<> per_pixel_value(rows, cols);
  cpp11::writable::doubles_matrix<> per_pixel_combined(rows, cols);
  cpp11::writable::doubles_matrix<> per_pixel_value_scaled(rows, cols);

  std::vector<std::vector<double>> centers_vals_vec(ncenters, std::vector<double>(bands));
  for (int i = 0; i < ncenters; i++) {
    for (int nval = 0; nval < bands; nval++) {
      centers_vals_vec[i][nval] = centers_vals(i, nval);
    }
  }

  std::vector<double> pixel_values;
  pixel_values.reserve(bands);

  std::vector<double> max_value_dist;
  if (adaptive_compactness) {
    max_value_dist = sc_compute_max_value_dist(centers_vals_vec, centers_xy, vals,
                                               rows, cols, bands, step, dist_name, dist_fun);
  }

  for (int i = 0; i < cols; i++) {
    for (int j = 0; j < rows; j++) {
      int cid = clusters(j, i);
      if (cid < 0 || cid >= ncenters) {
        per_pixel_spatial(j, i) = NA_REAL;
        per_pixel_value(j, i) = NA_REAL;
        per_pixel_combined(j, i) = NA_REAL;
        per_pixel_value_scaled(j, i) = NA_REAL;
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
        per_pixel_spatial(j, i) = NA_REAL;
        per_pixel_value(j, i) = NA_REAL;
        per_pixel_combined(j, i) = NA_REAL;
        per_pixel_value_scaled(j, i) = NA_REAL;
        continue;
      }

      /* value distance */
      double value_dist = get_vals_dist(centers_vals_vec[cid], pixel_values, dist_name, dist_fun);

      /* spatial distance */
      double center_x = centers_xy(cid, 0);
      double center_y = centers_xy(cid, 1);
      double y_dist = center_y - j;
      double x_dist = center_x - i;
      double spatial_dist = sqrt((y_dist * y_dist) + (x_dist * x_dist));

      double denom = compactness;
      if (adaptive_compactness) {
        denom = max_value_dist[cid];
      }

      /* combined distance */
      double combined_dist = NA_REAL;
      if (denom != 0.0 && step != 0) {
        double dist1 = value_dist / denom;
        double dist2 = spatial_dist / step;
        combined_dist = sqrt((dist1 * dist1) + (dist2 * dist2));
      }

      per_pixel_spatial(j, i) = spatial_dist;
      per_pixel_value(j, i) = value_dist;
      per_pixel_combined(j, i) = combined_dist;
      per_pixel_value_scaled(j, i) = (denom != 0.0) ? value_dist / denom : NA_REAL;
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
