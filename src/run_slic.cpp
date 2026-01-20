#include "slic.h"

[[cpp11::register]]
list run_slic(cpp11::integers mat, cpp11::doubles_matrix<> vals, int step, double compactness, bool clean, bool centers,
              std::string dist_name, cpp11::function dist_fun,
              cpp11::function avg_fun_fun, std::string avg_fun_name, int iter, int minarea,
              cpp11::integers_matrix<> input_centers, int verbose, bool diagnostics) {
  /* Show the step value */
  if (verbose > 0) Rprintf("Step: %u\n", step);

  /* Create an object of the slic class. */
  Slic slic;

  /* Create supercells */
  slic.generate_superpixels(mat, vals, step, compactness, dist_name, dist_fun, avg_fun_fun, avg_fun_name, iter, input_centers, verbose, diagnostics);

  /* Enforce connectivity if wanted */
  if (clean){
    slic.create_connectivity(vals, avg_fun_fun, avg_fun_name, minarea, verbose);
  }

  /* Write the result as a list */
  writable::list result(4);
  result.at(0) = slic.return_clusters();
  if (centers) {
    result.at(1) = slic.return_centers();
    result.at(2) = slic.return_centers_vals();
  } else {
    result.at(1) = R_NilValue;
    result.at(2) = R_NilValue;
  }
  if (diagnostics) {
    result.at(3) = slic.build_diagnostics(vals, dist_name, dist_fun);
  } else {
    result.at(3) = R_NilValue;
  }

  /* Return the result */
  return result;
}
