#pragma once
#include "cpp11.hpp"
#include <iostream>
#include <vector>
#include <float.h>

using namespace std;
using namespace cpp11;
namespace writable = cpp11::writable;

#define NR_ITERATIONS 10

class Slic {
  private:
    /* The cluster assignments and distance values for each pixel. */
    vector<vector<int> > clusters;
    vector<vector<double> > distances;

    /* The values of the centers. */
    vector<vector<double> > centers;
    vector<vector<double> > centers_vals;

    /* The number of occurrences of each center. */
    vector<int> center_counts;

    /* Initialize the new cluster matrix. */
    vector<vector<int> > new_clusters;

    /* The step size per cluster, and the color (nc) and distance (ns)
     * parameters. */
    int step, nc, ns;

    double euclidean(vector<double> values1, vector<double> values2);
    double manhattan(vector<double> values1, vector<double> values2);
    double jensen_shannon(vector<double> values1, vector<double> values2);

    double get_vals_dist(vector<double> values1, vector<double> values2, std::string type);

    /* Compute the distance between a center and an individual pixel. */
    double compute_dist(int ci, int y, int x, vector<double> value, std::string type);

    /* Find the pixel with the lowest gradient in a 3x3 surrounding. */
    vector<int> find_local_minimum(doubles_matrix vals, integers_matrix mat, int y, int x, std::string type);

    /* Remove and initialize the 2d vectors. */
    void clear_data();
    void inits(integers_matrix m, doubles_matrix vals, std::string type);

  public:
    /* Class constructors and deconstructors. */
    Slic();
    ~Slic();

    /* Generate an over-segmentation for an image. */
    void generate_superpixels(integers_matrix m, doubles_matrix vals, double step, int nc, std::string type);
    /* Enforce connectivity for an image. */
    void create_connectivity(integers_matrix mat, doubles_matrix vals);

    writable::integers_matrix return_centers();
    writable::integers_matrix return_clusters();
};

[[cpp11::register]]
integers_matrix run_slic(integers_matrix m, doubles_matrix vals, double step, int nc, bool con, bool output_type, std::string type) {
  Slic slic;
  slic.generate_superpixels(m, vals, step, nc, type);
  if (con){
    slic.create_connectivity(m, vals);
  }
  if (output_type) {
    integers_matrix result = slic.return_clusters();
    return result;
  } else {
    integers_matrix result = slic.return_centers();
    return result;
  }
}
