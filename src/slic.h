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

typedef multimap<int, int> IntToIntMap;
typedef IntToIntMap::iterator mapIter;

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
    int step;
    double compactness;

    /* Create centers of initial supercells, adds their values, etc. */
    void create_centers(vector<int> mat_dims, doubles_matrix<> vals,
                        std::string& dist_name, cpp11::function dist_fun,
                        double step);

    /* Create centers of initial supercells based on the input provided by a user, adds their values, etc. */
    void create_centers2(vector<int> mat_dims, doubles_matrix<> vals,
                         std::string& dist_name, cpp11::function dist_fun,
                         integers_matrix<> input_centers);

    /* Compute the total (spatial and value) distance between a center and an individual pixel. */
    double compute_dist(int& ci, int& y, int& x, vector<double>& value,
                        std::string& dist_name, cpp11::function dist_fun);

    /* Find the pixel with the lowest gradient in a 3x3 surrounding. */
    vector<int> find_local_minimum(doubles_matrix<> vals, int& y, int& x,
                                   std::string& dist_name, cpp11::function dist_fun);

    /* Remove and initialize the 2d vectors. */
    void clear_data();

    /* Initialize the clusters, distances, and centers. */
    void inits(integers mat, doubles_matrix<> vals,
               std::string& type, cpp11::function dist_fun,
               integers_matrix<> input_centers);

  public:
    /* Class constructors and deconstructors. */
    Slic();
    ~Slic();

    /* Generate an over-segmentation for an image. */
    void generate_superpixels(integers mat, doubles_matrix<> vals, double step, double compactness,
                              std::string& dist_name, cpp11::function dist_fun,
                              cpp11::function avg_fun_fun, std::string& avg_fun_name, int iter,
                              integers_matrix<> input_centers, int verbose);

    /* Enforce connectivity for an image. */
    void create_connectivity(doubles_matrix<> vals, cpp11::function avg_fun_fun, std::string& avg_fun_name, int minarea, int verbose);

    /* Writes the results to a format usable by R/cpp11. */
    /* Returns final pixel assignments */
    writable::integers_matrix<> return_centers();
    /* Returns supercells centers*/
    writable::integers_matrix<> return_clusters();
    /* Returns supercells centers' values*/
    writable::doubles_matrix<> return_centers_vals();
};
