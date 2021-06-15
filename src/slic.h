#pragma once
#include "cpp11.hpp"
#include <iostream>
#include <vector>
#include <float.h>

using namespace std;
using namespace cpp11;
namespace writable = cpp11::writable;

// #define NR_ITERATIONS 10

class Slic {
  private:
    /* The cluster assignments and distance values for each pixel. */
    vector<vector<int> > clusters;
    vector<vector<double> > distances;

    /* Input matrix dimensions (rows, cols, val-cols) */
    vector<int> mat_dims;

    /* The values of the centers. */
    vector<vector<double> > centers;
    vector<vector<double> > centers_vals;

    /* The number of occurrences of each center. */
    vector<int> center_counts;

    /* Initialize the new cluster matrix. */
    vector<vector<int> > new_clusters;

    /* The step size per cluster, and the color (nc) and distance (ns)
     * parameters. */
    int step, ns;
    double nc;

    // double euclidean(vector<double>& values1, vector<double>& values2);
    // double manhattan(vector<double>& values1, vector<double>& values2);
    // double jensen_shannon(vector<double>& values1, vector<double>& values2);

    double get_vals_dist(vector<double>& values1, vector<double>& values2, std::string& type);

    /* Compute the distance between a center and an individual pixel. */
    double compute_dist(int& ci, int& y, int& x, vector<double>& value, std::string& type);

    /* Find the pixel with the lowest gradient in a 3x3 surrounding. */
    vector<int> find_local_minimum(doubles_matrix vals, int& y, int& x, std::string& type);

    /* Remove and initialize the 2d vectors. */
    void clear_data();
    void inits(integers mat, doubles_matrix vals, std::string& type);

  public:
    /* Class constructors and deconstructors. */
    Slic();
    ~Slic();

    /* Generate an over-segmentation for an image. */
    void generate_superpixels(integers mat, doubles_matrix vals, double step, double nc, std::string& type, int iter);
    /* Enforce connectivity for an image. */
    void create_connectivity(doubles_matrix vals);

    writable::doubles_matrix return_centers();
    writable::doubles_matrix return_centers_vals();
    writable::integers_matrix return_clusters();
};

[[cpp11::register]]
list run_slic(integers mat, doubles_matrix vals, int k, double nc, bool con, bool centers, std::string type, int iter) {

  const int superpixelsize = 0.5 + double(mat.at(0) * mat.at(1)) / double(k);
  const int step = sqrt(double(superpixelsize)) + 0.5;
  // cout << "superpixelsize" << superpixelsize << endl;
  // cout << "step" << step << endl;

  Slic slic;
  slic.generate_superpixels(mat, vals, step, nc, type, iter);
  if (con){
    slic.create_connectivity(vals);
  }
  writable::list result(3);
  result.at(0) = slic.return_clusters();
  if (centers) {
    result.at(1) = slic.return_centers();
    result.at(2) = slic.return_centers_vals();
  }
  return result;
}
