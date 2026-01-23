#pragma once
#include <array>
#include <cstdio>
#include <cmath>
#include <float.h>
#include <functional>
#include <map>
#include <numeric>
#include <string>
#include <vector>

class SlicCore {
 public:
  using DistFn = std::function<double(const std::vector<double>&, const std::vector<double>&)>;
  using AvgFn = std::function<double(const std::vector<double>&)>;

  SlicCore();
  ~SlicCore();

  void generate_superpixels(const std::vector<int>& mat_dims, const std::vector<double>& vals,
                            int step, double compactness, DistFn dist_fn, AvgFn avg_fn,
                            const std::string& avg_fun_name, int iter,
                            const std::vector<std::array<int, 2>>& input_centers,
                            int verbose, bool iter_diagnostics);

  void create_connectivity(const std::vector<double>& vals, AvgFn avg_fn,
                           const std::string& avg_fun_name, int minarea, int verbose);

  const std::vector<std::vector<int>>& clusters_ref() const { return clusters; }
  const std::vector<std::vector<double>>& centers_ref() const { return centers; }
  const std::vector<std::vector<double>>& centers_vals_ref() const { return centers_vals; }
  const std::vector<int>& mat_dims_ref() const { return mat_dims; }

  const std::vector<double>& iter_mean_distance_ref() const { return iter_mean_distance; }

  int step_value() const { return step; }
  double compactness_value() const { return compactness; }

 private:
  using IntToIntMap = std::multimap<int, int>;
  using mapIter = IntToIntMap::iterator;

  std::vector<std::vector<int>> clusters;
  std::vector<std::vector<double>> distances;
  std::vector<int> mat_dims;
  std::vector<std::vector<double>> centers;
  std::vector<std::vector<double>> centers_vals;
  std::vector<int> center_counts;
  std::vector<std::vector<int>> new_clusters;

  bool iter_diagnostics_enabled = false;
  std::vector<double> iter_mean_distance;

  int step = 0;
  double compactness = 0.0;

  const std::vector<double>* vals_ptr = nullptr;
  int bands = 0;
  DistFn dist_fn;
  AvgFn avg_fn;
  std::string avg_fun_name;

  void create_centers(const std::vector<int>& mat_dims, const std::vector<double>& vals, int step);
  void create_centers2(const std::vector<int>& mat_dims, const std::vector<double>& vals,
                       const std::vector<std::array<int, 2>>& input_centers);
  void clear_data();
  void inits(const std::vector<int>& mat_dims, const std::vector<double>& vals,
             const std::vector<std::array<int, 2>>& input_centers);

  double compute_dist(int ci, int y, int x, const std::vector<double>& values) const;
  std::vector<double> find_local_minimum(const std::vector<double>& vals, int y, int x);
  double value_at(int ncell, int nval) const;
};
