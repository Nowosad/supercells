#include "slic_core.h"
#include "distances.h"
#include "cpp11.hpp"
#include "cpp11/list.hpp"
#include "cpp11/matrix.hpp"

cpp11::writable::integers_matrix<> return_clusters(const SlicCore& slic);
cpp11::writable::doubles_matrix<> return_centers(const SlicCore& slic);
cpp11::writable::doubles_matrix<> return_centers_vals(const SlicCore& slic);
[[cpp11::register]]
cpp11::list run_slic(cpp11::integers mat, cpp11::doubles_matrix<> vals, int step, double compactness,
                     bool clean, bool centers, std::string dist_name, cpp11::function dist_fun,
                     cpp11::function avg_fun_fun, std::string avg_fun_name, int iter, int minarea,
                     cpp11::integers_matrix<> input_centers, int verbose, bool iter_diagnostics) {
  if (verbose > 0) Rprintf("Step: %u\n", step);

  int ncell = vals.nrow();
  int bands = vals.ncol();
  std::vector<double> vals_vec(ncell * bands);
  for (int i = 0; i < ncell; i++) {
    for (int j = 0; j < bands; j++) {
      vals_vec[i * bands + j] = vals(i, j);
    }
  }

  std::vector<int> mat_dims = {mat.at(0), mat.at(1), bands};

  std::vector<std::array<int, 2>> centers_vec;
  centers_vec.reserve(input_centers.nrow());
  for (int i = 0; i < input_centers.nrow(); i++) {
    centers_vec.push_back({input_centers(i, 0), input_centers(i, 1)});
  }

  auto dist_cb = [&](const std::vector<double>& values1, const std::vector<double>& values2) -> double {
    return get_vals_dist(values1, values2, dist_name, dist_fun);
  };

  SlicCore::AvgFn avg_cb = nullptr;
  if (avg_fun_name.empty()) {
    avg_cb = [&](const std::vector<double>& values) -> double {
      return cpp11::as_cpp<double>(avg_fun_fun(values));
    };
  }

  SlicCore slic;
  slic.generate_superpixels(mat_dims, vals_vec, step, compactness, dist_cb, avg_cb,
                            avg_fun_name, iter, centers_vec, verbose, iter_diagnostics);

  if (clean) {
    slic.create_connectivity(vals_vec, avg_cb, avg_fun_name, minarea, verbose);
  }

  cpp11::writable::list result(4);
  result.at(0) = return_clusters(slic);
  if (centers) {
    result.at(1) = return_centers(slic);
    result.at(2) = return_centers_vals(slic);
  } else {
    result.at(1) = R_NilValue;
    result.at(2) = R_NilValue;
  }
  if (iter_diagnostics) {
    cpp11::writable::list iter(1);
    iter.names() = {"mean_distance"};
    cpp11::writable::doubles iter_mean(slic.iter_mean_distance_ref().size());
    for (size_t i = 0; i < slic.iter_mean_distance_ref().size(); i++) {
      iter_mean[i] = slic.iter_mean_distance_ref()[i];
    }
    iter.at(0) = iter_mean;
    result.at(3) = iter;
  } else {
    result.at(3) = R_NilValue;
  }
  return result;
}
