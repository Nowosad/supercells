#include "calc_stats.h"

double median(vector<double>& v){
  size_t n = v.size() / 2;
  std::nth_element(v.begin(), v.begin()+n, v.end());
  return v[n];
}

double mean(vector<double>& v){
  double sum = std::accumulate(v.begin(), v.end(), 0.0);
  double mean = sum / v.size();
  return mean;
}
