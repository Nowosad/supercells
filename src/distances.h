#include "cpp11.hpp"
#include <iostream>
#include <vector>
#include <float.h>

using namespace std;

/* Calculate the values distance between a center and an individual pixel. */
double get_vals_dist(vector<double>& values1, vector<double>& values2,
                           std::string& dist_name, cpp11::function dist_fun);
double euclidean(std::vector<double>& values1, std::vector<double>& values2);
double manhattan(vector<double>& values1, vector<double>& values2);
double jensen_shannon(vector<double>& values1, vector<double>& values2);
double custom_log2(const double& x);
double dtw3(std::vector<double>& values1, std::vector<double>& values2);
double dtw2d(std::vector<double>& values1, std::vector<double>& values2);
double custom_distance(vector<double>& values1, vector<double>& values2, std::string& dist_name);
