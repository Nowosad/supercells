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
    /* The number of occurrences of each center. */
    vector<int> center_counts;

    /* Initialize the new cluster matrix. */
    vector<vector<int> > new_clusters;

    /* The step size per cluster, and the color (nc) and distance (ns)
     * parameters. */
    int step, nc, ns;

    /* Compute the distance between a center and an individual pixel. */
    double compute_dist(int ci, int y, int x, int value);

    /* Find the pixel with the lowest gradient in a 3x3 surrounding. */
    vector<int> find_local_minimum(integers_matrix m, int y, int x);

    /* Remove and initialize the 2d vectors. */
    void clear_data();
    void inits(integers_matrix m);

  public:
    /* Class constructors and deconstructors. */
    Slic();
    ~Slic();

    /* Generate an over-segmentation for an image. */
    writable::integers_matrix generate_superpixels(integers_matrix m, int step, int nc);
    /* Enforce connectivity for an image. */
    void create_connectivity(integers_matrix mat);

    writable::integers_matrix return_centers();
    writable::integers_matrix return_clusters();
};

[[cpp11::register]]
integers_matrix run_slic(integers_matrix m, int step, int nc, int con, int type) {
  Slic slic;
  integers_matrix r = slic.generate_superpixels(m, step, nc);
  if (con){
    slic.create_connectivity(m);
  }
  if (type) {
    integers_matrix result = slic.return_clusters();
    return result;
  } else {
    integers_matrix result = slic.return_centers();
    return result;
  }
}
