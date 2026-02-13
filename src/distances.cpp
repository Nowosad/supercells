#include "distances.h"
// using namespace cpp11::literals; // so we can use ""_nm syntax

double get_vals_dist(const std::vector<double>& values1, const std::vector<double>& values2,
                     const std::string& dist_name, cpp11::function dist_fun){
  if (dist_name == "euclidean"){
    return euclidean(values1, values2);
  } else if (dist_name == "jsd"){
    return jensen_shannon(values1, values2);
  } else if (dist_name == "dtw"){
    return dtw1d(values1, values2);
  } else if (dist_name == "dtw2d"){
    return dtw2d(values1, values2);
  } else if (dist_name != ""){
    return custom_distance(values1, values2, dist_name);
  } else {
    return cpp11::as_cpp<double>(dist_fun(values1, values2));
  }
}

double euclidean(const std::vector<double>& values1, const std::vector<double>& values2){
  int len1 = values1.size();
  double dist = 0.0;
  double diff = 0.0;

  for (int i = 0; i < len1; i++){
    diff = values1[i] - values2[i];
    dist += diff * diff;
  }
  return sqrt(dist);
}

double manhattan(const std::vector<double>& values1, const std::vector<double>& values2){
  int len1 = values1.size();
  double dist = 0.0;
  double diff = 0.0;

  for (int i = 0; i < len1; i++){
    diff = fabs(values1[i] - values2[i]);
    dist += diff;
  }
  return dist;
}

double dtw1d(const std::vector<double>& values1, const std::vector<double>& values2){
  int len1 = values1.size();
  int len2 = values2.size();
  if (len1 == 0 || len2 == 0) {
    return NAN;
  }

  std::vector<double> prev(len2);
  std::vector<double> curr(len2);

  prev[0] = std::fabs(values1[0] - values2[0]);
  for (int j = 1; j < len2; j++) {
    prev[j] = prev[j - 1] + std::fabs(values1[0] - values2[j]);
  }

  for (int i = 1; i < len1; i++) {
    curr[0] = prev[0] + std::fabs(values1[i] - values2[0]);
    for (int j = 1; j < len2; j++) {
      double cost = std::fabs(values1[i] - values2[j]);
      double best = std::fmin(std::fmin(prev[j], curr[j - 1]), prev[j - 1]);
      curr[j] = cost + best;
    }
    prev.swap(curr);
  }

  return prev[len2 - 1];
}

double dtw2d(const std::vector<double>& values1, const std::vector<double>& values2){
  int len1 = values1.size() / 2;
  int len2 = values2.size() / 2;
  if (len1 == 0 || len2 == 0) {
    return NAN;
  }

  std::vector<double> prev(len2);
  std::vector<double> curr(len2);

  auto cost = [&](int i, int j) {
    double dx = values1[i + len1] - values2[j + len2];
    double dy = values1[i] - values2[j];
    return std::sqrt(dx * dx + dy * dy);
  };

  prev[0] = cost(0, 0);
  for (int j = 1; j < len2; j++) {
    prev[j] = prev[j - 1] + cost(0, j);
  }

  for (int i = 1; i < len1; i++) {
    curr[0] = prev[0] + cost(i, 0);
    for (int j = 1; j < len2; j++) {
      double best = std::fmin(std::fmin(prev[j], curr[j - 1]), prev[j - 1]);
      curr[j] = cost(i, j) + best;
    }
    prev.swap(curr);
  }

  return prev[len2 - 1];
}


double jensen_shannon(const std::vector<double>& values1, const std::vector<double>& values2){
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

double custom_log2(double x){
  if (x == 0.0){
    return NAN;
  } else {
    return log(x)/log(2.0);
  }
}

double custom_distance(const std::vector<double>& values1, const std::vector<double>& values2, const std::string& dist_name){
  static cpp11::function single_distance = cpp11::package("philentropy")["dist_one_one"];
  double p = NA_REAL;
  bool testNA = false;
  std::string unit = "log2";
  return cpp11::as_cpp<double>(single_distance(values1, values2, dist_name, p, testNA, unit));
}
