#include "distances.h"
#include "cpp11.hpp"
#include "cpp11/doubles.hpp"
#include <vector>

// Expose distance calculation to R for reuse in merging utilities.
[[cpp11::register]]
double sc_dist_vec_cpp(cpp11::doubles a,
                       cpp11::doubles b,
                       std::string dist_name,
                       cpp11::function dist_fun) {
  if (a.size() != b.size()) {
    cpp11::stop("Input vectors must have the same length");
  }
  std::vector<double> va(a.size());
  std::vector<double> vb(b.size());
  for (int i = 0; i < a.size(); i++) {
    va[i] = a[i];
    vb[i] = b[i];
  }
  return get_vals_dist(va, vb, dist_name, dist_fun);
}
