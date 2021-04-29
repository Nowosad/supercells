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

    /* The LAB and xy values of the centers. */
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


    /* Remove and initialize the 2d vectors. */
    void clear_data();
    void inits(integers_matrix m);

  public:
    /* Class constructors and deconstructors. */
    Slic();
    ~Slic();

    /* Find the pixel with the lowest gradient in a 3x3 surrounding. */
    vector<int> find_local_minimum(integers_matrix m, int y, int x);

    /* Generate an over-segmentation for an image. */
    writable::integers_matrix generate_superpixels(integers_matrix m, int step, int nc);
    /* Enforce connectivity for an image. */
    void create_connectivity(integers_matrix m);

    /* Draw functions. Resp. displayal of the centers and the contours. */
    // void display_center_grid(IplImage *image, CvScalar colour);
    // void display_contours(IplImage *image, CvScalar colour);
    // void colour_with_cluster_means(IplImage *image);
    writable::integers_matrix colour_with_cluster_means(integers_matrix m);
};

[[cpp11::register]]
integers_matrix foo(integers_matrix m, int step, int nc) {
  // integers_matrix r(20, 3);
  Slic slic;
  integers_matrix r = slic.generate_superpixels(m, step, nc);
  return r;
}

[[cpp11::register]]
integers_matrix foo2(integers_matrix m, int x, int y) {
  // integers_matrix r(20, 3);
  Slic slic;
  vector<int> r0 = slic.find_local_minimum(m, x, y);

  cout << "r0" << r0[0] << endl;

  writable::integers_matrix result(1, 2);
  result(0, 0) = r0[0];
  result(0, 1) = r0[1];

  return result;
}

[[cpp11::register]]
integers_matrix foo3(integers_matrix m, int step, int nc) {
  // integers_matrix r(20, 3);
  Slic slic;
  integers_matrix r = slic.generate_superpixels(m, step, nc);
  // slic.create_connectivity(m);
  integers_matrix r2 = slic.colour_with_cluster_means(m);
  return r2;
}
