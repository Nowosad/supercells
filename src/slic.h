#pragma once
#include "cpp11/matrix.hpp"
#include "cpp11/function.hpp"
#include "cpp11/integers.hpp"
#include "cpp11/list.hpp"
#include <iostream>
#include <vector>
#include <map>
#include <numeric>
#include <float.h>
#include "distances.h"
#include "calc_stats.h"

using namespace std;
using namespace cpp11;
// namespace writable = cpp11::writable;

// #define NR_ITERATIONS 10

class Slic {
  private:
    /* The cluster assignments and distance values for each pixel. */
    vector<vector<int> > clusters;
    vector<vector<double> > distances;

    /* Input matrix dimensions (rows, cols, val-cols) */
    vector<int> mat_dims;

    /* The values of the centers. */
    vector<vector<int> > centers;
    vector<vector<double> > centers_vals;

    /* The number of occurrences of each center. */
    vector<int> center_counts;

    /* Initialize the new cluster matrix. */
    vector<vector<int> > new_clusters;

    /* The step size per cluster, and the color (nc) and distance (ns)
     * parameters. */
    int step, ns;
    double nc;

    void create_centers(vector<int> mat_dims, doubles_matrix<> vals,
                        std::string& type, cpp11::function type_fun, double step);

    void create_centers2(vector<int> mat_dims,
                          doubles_matrix<> vals, std::string& type,
                          cpp11::function type_fun,
                          integers_matrix<> input_centers);

    double get_vals_dist(vector<double>& values1, vector<double>& values2,
                         std::string& type, cpp11::function type_fun);

    /* Compute the distance between a center and an individual pixel. */
    double compute_dist(int& ci, int& y, int& x, vector<double>& value,
                        std::string& type, cpp11::function type_fun);

    /* Find the pixel with the lowest gradient in a 3x3 surrounding. */
    vector<int> find_local_minimum(doubles_matrix<> vals, int& y, int& x,
                                   std::string& type, cpp11::function type_fun);

    /* Remove and initialize the 2d vectors. */
    void clear_data();
    void inits(integers mat, doubles_matrix<> vals,
               std::string& type, cpp11::function type_fun, integers_matrix<> input_centers);

  public:
    /* Class constructors and deconstructors. */
    Slic();
    ~Slic();

    /* Generate an over-segmentation for an image. */
    void generate_superpixels(integers mat, doubles_matrix<> vals, double step, double nc,
                              std::string& type, cpp11::function type_fun,
                              cpp11::function avg_fun_fun, std::string& avg_fun_name, int iter,
                              integers_matrix<> input_centers, int verbose);
    /* Enforce connectivity for an image. */
    void create_connectivity(doubles_matrix<> vals, cpp11::function avg_fun_fun, std::string& avg_fun_name, int lims, int verbose);

    writable::integers_matrix<> return_centers();
    writable::doubles_matrix<> return_centers_vals();
    writable::integers_matrix<> return_clusters();
};

[[cpp11::register]]
list run_slic(cpp11::integers mat, cpp11::doubles_matrix<> vals, int step, double nc, bool con, bool centers,
              std::string type, cpp11::function type_fun,
              cpp11::function avg_fun_fun, std::string avg_fun_name, int iter, int lims,
              cpp11::integers_matrix<> input_centers, int verbose) {

  // cout << "superpixelsize" << superpixelsize << endl;
  if (verbose > 0) Rprintf("Step: %u\n", step);

  // Rprintf("Vf: %f\n", vals(0, 0));

  Slic slic;
  slic.generate_superpixels(mat, vals, step, nc, type, type_fun, avg_fun_fun, avg_fun_name, iter, input_centers, verbose);
  if (con){
    slic.create_connectivity(vals, avg_fun_fun, avg_fun_name, lims, verbose);
  }
  writable::list result(3);
  result.at(0) = slic.return_clusters();
  if (centers) {
    result.at(1) = slic.return_centers();
    result.at(2) = slic.return_centers_vals();
  }
  return result;
}
