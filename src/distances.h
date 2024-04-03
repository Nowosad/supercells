#include "cpp11.hpp"
#include <iostream>
#include <vector>
#include <float.h>

using namespace std;

/* Calculate the values distance between a center and an individual pixel. */

double get_vals_dist(const vector<double>& values1, const vector<double>& values2,
                          const std::string dist_name, cpp11::function dist_fun);

[[cpp11::register]]
double get_vals_dist_cpp11(cpp11::doubles values1, cpp11::doubles values2,
                           const std::string dist_name, cpp11::function dist_fun);

double euclidean(const std::vector<double>& values1, const std::vector<double>& values2);
double manhattan(const vector<double>& values1, const vector<double>& values2);
double jensen_shannon(const vector<double>& values1, const vector<double>& values2);
double custom_log2(const double& x);
double dtw3(const std::vector<double>& values1, const std::vector<double>& values2);
double dtw2d(const std::vector<double>& values1, const std::vector<double>& values2);
double custom_distance(const vector<double>& values1, const vector<double>& values2, const std::string dist_name);
