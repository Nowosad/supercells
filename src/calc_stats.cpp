#include "calc_stats.h"

double median(std::vector<double>& v){
  size_t n = v.size() / 2;
  std::nth_element(v.begin(), v.begin()+n, v.end());
  double upper = v[n];
  if (v.size() % 2 != 0) return upper;
  double lower = *std::max_element(v.begin(), v.begin()+n);
  return lower + (upper - lower) / 2.0;
}

double mean(const std::vector<double>& v){
  double sum = std::accumulate(v.begin(), v.end(), 0.0);
  return sum / v.size();
}
