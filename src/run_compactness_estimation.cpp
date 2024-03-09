#include "slic.h"

[[cpp11::register]]
doubles run_compactness_estimation(integers mat, doubles_matrix<> vals, double step,
                                                     std::string dist_name, cpp11::function dist_fun,
                                                     integers_matrix<> input_centers){
  /* Create an object of the slic class. */
  Slic slic;

  doubles result = slic.estimate_compactness(mat, vals, step, dist_name, dist_fun, input_centers);
  // cout << result[0] << endl;
  return result;
}
