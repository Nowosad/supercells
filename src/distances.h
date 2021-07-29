#include "cpp11.hpp"
#include <iostream>
#include <vector>
#include <float.h>

using namespace std;

double euclidean(std::vector<double>& values1, std::vector<double>& values2);
double manhattan(vector<double>& values1, vector<double>& values2);
double jensen_shannon(vector<double>& values1, vector<double>& values2);
double custom_log2(const double& x);
double dtw3(std::vector<double>& values1, std::vector<double>& values2);
double custom_distance(vector<double>& values1, vector<double>& values2, std::string& type);
