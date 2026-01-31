#include "metrics_helpers.h"
#include "cpp11.hpp"
#include "cpp11/matrix.hpp"
#include <vector>

// Compute per-center mean value distance in a local 2*step window.
// Role: support local compactness estimation for sc_tune_compactness().
// Inputs: center coordinates/values, full raster values, window size, and distance function.
// Output: vector of mean distances per center (NA if no valid pixels).
[[cpp11::register]]
cpp11::writable::doubles sc_metrics_local_mean_cpp(cpp11::doubles_matrix<> centers_xy,
                                                   cpp11::doubles_matrix<> centers_vals,
                                                   cpp11::doubles_matrix<> vals,
                                                   int rows, int cols, int step,
                                                   std::string dist_name, cpp11::function dist_fun) {
  int bands = vals.ncol();
  int ncenters = centers_vals.nrow();

  std::vector<std::vector<double>> centers_vals_vec(ncenters, std::vector<double>(bands));
  for (int i = 0; i < ncenters; i++) {
    for (int nval = 0; nval < bands; nval++) {
      centers_vals_vec[i][nval] = centers_vals(i, nval);
    }
  }

  cpp11::writable::doubles result(ncenters);
  std::vector<double> pixel_values;
  pixel_values.reserve(bands);

  for (int i = 0; i < ncenters; i++) {
    double sum_dist = 0.0;
    int count = 0;
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
        sum_dist += d;
        count += 1;
      }
    }
    if (count > 0) {
      result[i] = sum_dist / count;
    } else {
      result[i] = NA_REAL;
    }
  }

  return result;
}
