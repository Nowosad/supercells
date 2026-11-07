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

  std::vector<std::vector<double>> center_values = sc_copy_center_values(centers_vals);

  cpp11::writable::doubles result(ncenters);
  std::vector<double> pixel_values;
  pixel_values.reserve(bands);

  for (int center = 0; center < ncenters; center++) {
    double distance_sum = 0.0;
    int pixel_count = 0;
    int center_col = static_cast<int>(std::round(centers_xy(center, 0)));
    int center_row = static_cast<int>(std::round(centers_xy(center, 1)));
    for (int col = center_col - step; col < center_col + step; col++) {
      for (int row = center_row - step; row < center_row + step; row++) {
        if (col < 0 || col >= cols || row < 0 || row >= rows) {
          continue;
        }
        int ncell = col + (row * cols);
        if (!sc_read_pixel_values(vals, ncell, pixel_values)) {
          continue;
        }
        double distance = get_vals_dist(
          center_values[center], pixel_values, dist_name, dist_fun
        );
        distance_sum += distance;
        pixel_count += 1;
      }
    }
    if (pixel_count > 0) {
      result[center] = distance_sum / pixel_count;
    } else {
      result[center] = NA_REAL;
    }
  }

  return result;
}
