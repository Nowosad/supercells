#pragma once

#include "distances.h"
#include "cpp11.hpp"
#include "cpp11/matrix.hpp"
#include <cmath>
#include <vector>

inline std::vector<std::vector<double>> sc_copy_center_values(
  const cpp11::doubles_matrix<>& centers_vals
) {
  int ncenters = centers_vals.nrow();
  int bands = centers_vals.ncol();
  std::vector<std::vector<double>> result(ncenters, std::vector<double>(bands));
  for (int center = 0; center < ncenters; center++) {
    for (int band = 0; band < bands; band++) {
      result[center][band] = centers_vals(center, band);
    }
  }
  return result;
}

inline bool sc_read_pixel_values(const cpp11::doubles_matrix<>& vals, int ncell,
                                 std::vector<double>& result) {
  result.clear();
  for (int band = 0; band < vals.ncol(); band++) {
    double value = vals(ncell, band);
    if (std::isnan(value)) {
      return false;
    }
    result.push_back(value);
  }
  return true;
}

// Compute per-center max value distance within a 2*step window.
// Role: provide local-window adaptive compactness denominators.
// Inputs: per-center values, center coordinates, full raster values, window size, distance function.
// Output: vector of max distances per center (min 1.0, NA pixels ignored).
inline std::vector<double> sc_compute_max_value_dist(
  const std::vector<std::vector<double>>& center_values,
  const cpp11::doubles_matrix<>& centers_xy,
  const cpp11::doubles_matrix<>& vals,
  int rows, int cols, int step,
  const std::string& dist_name,
  const cpp11::function& dist_fun
) {
  int ncenters = static_cast<int>(center_values.size());
  std::vector<double> max_value_dist(ncenters, 1.0);
  std::vector<double> pixel_values;
  pixel_values.reserve(vals.ncol());

  for (int center = 0; center < ncenters; center++) {
    double max_dist = 0.0;
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
        if (distance > max_dist) {
          max_dist = distance;
        }
      }
    }
    if (max_dist <= 0.0) {
      max_dist = 1.0;
    }
    max_value_dist[center] = max_dist;
  }
  return max_value_dist;
}

struct ScPixelMetrics {
  double value;
  double spatial;
  double combined;
  double value_scaled;
};

struct ScMetricData {
  std::vector<std::vector<double>> center_values;
  std::vector<double> denominators;
};

struct ScMetricSums {
  std::vector<double> value;
  std::vector<double> spatial;
  std::vector<double> combined;
  std::vector<int> count;
};

inline ScMetricData sc_prepare_metric_data(
  const cpp11::doubles_matrix<>& centers_xy,
  const cpp11::doubles_matrix<>& centers_vals,
  const cpp11::doubles_matrix<>& vals,
  int rows, int cols, int step, double compactness,
  bool adaptive_compactness, const std::string& dist_name,
  const cpp11::function& dist_fun
) {
  ScMetricData data;
  data.center_values = sc_copy_center_values(centers_vals);
  data.denominators.assign(data.center_values.size(), compactness);
  if (adaptive_compactness) {
    data.denominators = sc_compute_max_value_dist(
      data.center_values, centers_xy, vals, rows, cols, step, dist_name, dist_fun
    );
  }
  return data;
}

inline bool sc_compute_pixel_metrics(
  int cid, int col, int row, int cols, int step,
  const cpp11::doubles_matrix<>& centers_xy,
  const cpp11::doubles_matrix<>& vals,
  const ScMetricData& data,
  const std::string& dist_name,
  const cpp11::function& dist_fun,
  std::vector<double>& pixel_values,
  ScPixelMetrics& result
) {
  if (cid < 0 || cid >= static_cast<int>(data.center_values.size())) {
    return false;
  }

  int ncell = col + (row * cols);
  if (!sc_read_pixel_values(vals, ncell, pixel_values)) {
    return false;
  }

  result.value = get_vals_dist(data.center_values[cid], pixel_values, dist_name, dist_fun);

  double x_dist = centers_xy(cid, 0) - col;
  double y_dist = centers_xy(cid, 1) - row;
  result.spatial = std::sqrt((x_dist * x_dist) + (y_dist * y_dist));

  double denominator = data.denominators[cid];
  result.combined = NA_REAL;
  if (denominator != 0.0 && step != 0) {
    double value_scaled = result.value / denominator;
    double spatial_scaled = result.spatial / step;
    result.combined = std::sqrt(
      (value_scaled * value_scaled) + (spatial_scaled * spatial_scaled)
    );
  }
  result.value_scaled = (denominator != 0.0) ? result.value / denominator : NA_REAL;
  return true;
}

inline ScMetricSums sc_accumulate_metrics(
  const cpp11::integers_matrix<>& clusters,
  const cpp11::doubles_matrix<>& centers_xy,
  const cpp11::doubles_matrix<>& vals,
  const ScMetricData& data,
  int step, const std::string& dist_name,
  const cpp11::function& dist_fun
) {
  int ncenters = static_cast<int>(data.center_values.size());
  ScMetricSums sums{
    std::vector<double>(ncenters, 0.0),
    std::vector<double>(ncenters, 0.0),
    std::vector<double>(ncenters, 0.0),
    std::vector<int>(ncenters, 0)
  };
  std::vector<double> pixel_values;
  pixel_values.reserve(vals.ncol());

  for (int col = 0; col < clusters.ncol(); col++) {
    for (int row = 0; row < clusters.nrow(); row++) {
      int cid = clusters(row, col);
      ScPixelMetrics metrics;
      bool valid = sc_compute_pixel_metrics(
        cid, col, row, clusters.ncol(), step, centers_xy, vals, data,
        dist_name, dist_fun, pixel_values, metrics
      );
      if (!valid) {
        continue;
      }
      sums.value[cid] += metrics.value;
      sums.spatial[cid] += metrics.spatial;
      if (!std::isnan(metrics.combined)) {
        sums.combined[cid] += metrics.combined;
      }
      sums.count[cid] += 1;
    }
  }
  return sums;
}
