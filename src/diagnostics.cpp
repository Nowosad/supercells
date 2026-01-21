#include "slic_core.h"
#include "cpp11.hpp"
#include "cpp11/list.hpp"
#include "cpp11/matrix.hpp"
#include "cpp11/doubles.hpp"
#include <algorithm>
#include <cmath>
#include <numeric>

static double mean_vec(const std::vector<double>& v) {
  if (v.empty()) {
    return NA_REAL;
  }
  double sum = std::accumulate(v.begin(), v.end(), 0.0);
  return sum / v.size();
}

cpp11::writable::list build_iter_diagnostics(const SlicCore& slic) {
  cpp11::writable::list iteration(1);
  iteration.names() = {"mean_distance"};
  cpp11::writable::doubles iter_mean(slic.iter_mean_distance_ref().size());
  for (size_t i = 0; i < slic.iter_mean_distance_ref().size(); i++) {
    iter_mean[i] = slic.iter_mean_distance_ref()[i];
  }
  iteration.at(0) = iter_mean;

  return iteration;
}

cpp11::writable::list build_diagnostics(const SlicCore& slic, const std::vector<double>& vals,
                                        SlicCore::DistFn dist_fn) {
  const auto& mat_dims = slic.mat_dims_ref();
  const auto& clusters = slic.clusters_ref();
  const auto& centers = slic.centers_ref();
  const auto& centers_vals = slic.centers_vals_ref();

  int rows = mat_dims[0];
  int cols = mat_dims[1];
  int bands = mat_dims[2];
  int ncenters = centers.size();

  cpp11::writable::doubles_matrix<> per_pixel_spatial(rows, cols);
  cpp11::writable::doubles_matrix<> per_pixel_value(rows, cols);
  cpp11::writable::doubles_matrix<> per_pixel_combined(rows, cols);

  std::vector<std::vector<double> > cluster_spatial(ncenters);
  std::vector<std::vector<double> > cluster_value(ncenters);
  std::vector<double> all_spatial;
  std::vector<double> all_value;
  all_spatial.reserve(rows * cols);
  all_value.reserve(rows * cols);

  std::vector<double> pixel_values;
  pixel_values.reserve(bands);

  for (int i = 0; i < cols; i++) {
    for (int j = 0; j < rows; j++) {
      int cid = clusters[i][j];
      if (cid < 0) {
        per_pixel_spatial(j, i) = NA_REAL;
        per_pixel_value(j, i) = NA_REAL;
        per_pixel_combined(j, i) = NA_REAL;
        continue;
      }

      int ncell = i + (j * cols);
      pixel_values.clear();
      bool has_na = false;
      for (int nval = 0; nval < bands; nval++) {
        double val = vals[ncell * bands + nval];
        pixel_values.push_back(val);
        if (std::isnan(val)) {
          has_na = true;
        }
      }
      if (has_na) {
        per_pixel_spatial(j, i) = NA_REAL;
        per_pixel_value(j, i) = NA_REAL;
        per_pixel_combined(j, i) = NA_REAL;
        continue;
      }

      double value_dist = dist_fn(centers_vals[cid], pixel_values);
      int y_dist = centers[cid][0] - j;
      int x_dist = centers[cid][1] - i;
      double spatial_dist = sqrt((y_dist * y_dist) + (x_dist * x_dist));

      double combined_dist = NA_REAL;
      if (slic.compactness_value() != 0.0 && slic.step_value() != 0.0) {
        double dist1 = value_dist / slic.compactness_value();
        double dist2 = spatial_dist / slic.step_value();
        combined_dist = sqrt((dist1 * dist1) + (dist2 * dist2));
      }

      per_pixel_spatial(j, i) = spatial_dist;
      per_pixel_value(j, i) = value_dist;
      per_pixel_combined(j, i) = combined_dist;

      cluster_spatial[cid].push_back(spatial_dist);
      cluster_value[cid].push_back(value_dist);
      all_spatial.push_back(spatial_dist);
      all_value.push_back(value_dist);
    }
  }

  cpp11::writable::doubles mean_value(ncenters);
  cpp11::writable::doubles mean_spatial(ncenters);
  cpp11::writable::doubles ratio_mean(ncenters);
  for (int i = 0; i < ncenters; i++) {
    double mv = mean_vec(cluster_value[i]);
    double ms = mean_vec(cluster_spatial[i]);
    mean_value[i] = mv;
    mean_spatial[i] = ms;
    if (slic.compactness_value() != 0.0 && slic.step_value() != 0.0 && ms > 0.0) {
      ratio_mean[i] = (mv / slic.compactness_value()) / (ms / slic.step_value());
    } else {
      ratio_mean[i] = NA_REAL;
    }
  }

  std::vector<double> weighted_value;
  std::vector<double> weighted_spatial;
  if (slic.compactness_value() != 0.0) {
    weighted_value.reserve(all_value.size());
    for (size_t i = 0; i < all_value.size(); i++) {
      weighted_value.push_back(all_value[i] / slic.compactness_value());
    }
  }
  if (slic.step_value() != 0.0) {
    weighted_spatial.reserve(all_spatial.size());
    for (size_t i = 0; i < all_spatial.size(); i++) {
      weighted_spatial.push_back(all_spatial[i] / slic.step_value());
    }
  }

  double spatial_mean = mean_vec(all_spatial);
  double value_mean = mean_vec(all_value);
  double w_spatial_mean = mean_vec(weighted_spatial);
  double w_value_mean = mean_vec(weighted_value);

  cpp11::writable::list iteration = build_iter_diagnostics(slic);

  cpp11::writable::list per_cluster(3);
  per_cluster.names() = {"mean_value", "mean_spatial", "ratio_mean"};
  per_cluster.at(0) = mean_value;
  per_cluster.at(1) = mean_spatial;
  per_cluster.at(2) = ratio_mean;

  cpp11::writable::list per_pixel(3);
  per_pixel.names() = {"spatial", "value", "combined"};
  per_pixel.at(0) = per_pixel_spatial;
  per_pixel.at(1) = per_pixel_value;
  per_pixel.at(2) = per_pixel_combined;

  cpp11::writable::list scale(4);
  scale.names() = {"spatial_mean", "value_mean", "weighted_spatial_mean", "weighted_value_mean"};
  scale.at(0) = cpp11::as_sexp(spatial_mean);
  scale.at(1) = cpp11::as_sexp(value_mean);
  scale.at(2) = cpp11::as_sexp(w_spatial_mean);
  scale.at(3) = cpp11::as_sexp(w_value_mean);

  cpp11::writable::list result(4);
  result.names() = {"iteration", "per_cluster", "per_pixel", "scale"};
  result.at(0) = iteration;
  result.at(1) = per_cluster;
  result.at(2) = per_pixel;
  result.at(3) = scale;

  return result;
}
