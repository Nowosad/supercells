#include "distances.h"
#include "dtw/include/DTW.hpp"
// using namespace cpp11::literals; // so we can use ""_nm syntax

double get_vals_dist(vector<double>& values1, vector<double>& values2,
                           std::string& dist_name, cpp11::function dist_fun){
  if (dist_name == "euclidean"){
    return euclidean(values1, values2);
  } else if (dist_name == "jsd"){
    return jensen_shannon(values1, values2);
  } else if (dist_name == "dtw"){
    return dtw3(values1, values2);
  } else if (dist_name == "dtw2d"){
    return dtw2d(values1, values2);
  } else if (dist_name != ""){
    return custom_distance(values1, values2, dist_name);
  } else {
    return dist_fun(values1, values2);
  }
}

double euclidean(std::vector<double>& values1, std::vector<double>& values2){
  int len1 = values1.size();
  double dist = 0.0;
  double diff = 0.0;

  for (int i = 0; i < len1; i++){
    diff = values1[i] - values2[i];
    dist += diff * diff;
  }
  return sqrt(dist);
}

double manhattan(vector<double>& values1, vector<double>& values2){
  int len1 = values1.size();
  double dist = 0.0;
  double diff = 0.0;

  for (int i = 0; i < len1; i++){
    diff = fabs(values1[i] - values2[i]);
    dist += diff;
  }
  return dist;
}

double dtw3(std::vector<double>& values1, std::vector<double>& values2){
  double p = 2;  // the p-norm to use; 2.0 = euclidean, 1.0 = manhattan
  int len1 = values1.size();

  std::vector<std::vector<double> > a;
  std::vector<std::vector<double> > b;

  a.reserve(len1);
  b.reserve(len1);

  for(int i = 0; i < len1; i++){
    std::vector<double> pair1(2); std::vector<double> pair2(2);
    pair1[0] = i;
    pair1[1] = values1[i];
    pair2[0] = i;
    pair2[1] = values2[i];
    a.push_back(pair1); b.push_back(pair2);
  }

  double scost = DTW::dtw_distance_only(a, b, p);
  return(scost);
}

double dtw2d(std::vector<double>& values1, std::vector<double>& values2){
  double p = 2;  // the p-norm to use; 2.0 = euclidean, 1.0 = manhattan
  int len1 = values1.size() / 2;

  std::vector<std::vector<double> > a;
  std::vector<std::vector<double> > b;

  a.reserve(len1);
  b.reserve(len1);

  for(int i = 0; i < len1; i++){
    std::vector<double> pair1(2); std::vector<double> pair2(2);
    pair1[0] = values1[i + (len1 - 1)];
    pair1[1] = values1[i];
    pair2[0] = values2[i + (len1 - 1)];
    pair2[1] = values2[i];
    a.push_back(pair1); b.push_back(pair2);
  }

  double scost = DTW::dtw_distance_only(a, b, p);
  return(scost);
}


double jensen_shannon(vector<double>& values1, vector<double>& values2){
  int len1 = values1.size();
  double sum1       = 0.0;
  double sum2       = 0.0;
  double sum12      = 0.0;

  for(int i = 0; i < len1; i++){
    sum12 = values1[i] + values2[i];
    if (values1[i] == 0.0 || sum12 == 0.0) {
      sum1  += 0.0;
    } else {
      sum1  +=  values1[i] * custom_log2((2.0 * values1[i]) / sum12);
    }
    if (values2[i] == 0.0 || sum12 == 0.0) {
      sum2  += 0.0;
    } else {
      sum2  +=  values2[i] * custom_log2((2.0 * values2[i]) / sum12);
    }
  }
  return 0.5 * (sum1 + sum2);
}

double custom_log2(const double& x){
  if (x == 0.0){
    return NAN;
  } else {
    return log(x)/log(2.0);
  }
}

double custom_distance(vector<double>& values1, vector<double>& values2, std::string& dist_name){
  auto single_distance = cpp11::package("philentropy")["dist_one_one"];
  double p = NA_REAL;
  bool testNA = false;
  std::string unit = "log2";
  return single_distance(values1, values2, dist_name, p, testNA, unit);
}
