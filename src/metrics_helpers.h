#pragma once

#include "distances.h"
#include "cpp11.hpp"
#include "cpp11/matrix.hpp"
#include <cmath>
#include <vector>

// Compute per-center max value distance within a 2*step window.
// Role: provide SLIC0-style adaptive compactness denominators.
// Inputs: per-center values, center coordinates, full raster values, window size, distance function.
// Output: vector of max distances per center (min 1.0, NA pixels ignored).
inline std::vector<double> sc_compute_max_value_dist(
  const std::vector<std::vector<double>>& centers_vals_vec,
  const cpp11::doubles_matrix<>& centers_xy,
  const cpp11::doubles_matrix<>& vals,
  int rows, int cols, int bands, int step,
  const std::string& dist_name,
  const cpp11::function& dist_fun
) {
  int ncenters = static_cast<int>(centers_vals_vec.size());
  std::vector<double> max_value_dist(ncenters, 1.0);
  std::vector<double> pixel_values;
  pixel_values.reserve(bands);

  for (int i = 0; i < ncenters; i++) {
    double max_dist = 0.0;
    int center_x = static_cast<int>(std::round(centers_xy(i, 0)));
    int center_y = static_cast<int>(std::round(centers_xy(i, 1)));
    for (int m = center_x - step; m < center_x + step; m++) {
      for (int n = center_y - step; n < center_y + step; n++) {
        if (m < 0 || m >= cols || n < 0 || n >= rows) {
          continue;
        }
        int ncell = m + (n * cols);
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
        double d = get_vals_dist(centers_vals_vec[i], pixel_values, dist_name, dist_fun);
        if (d > max_dist) {
          max_dist = d;
        }
      }
    }
    if (max_dist <= 0.0) {
      max_dist = 1.0;
    }
    max_value_dist[i] = max_dist;
  }
  return max_value_dist;
}
