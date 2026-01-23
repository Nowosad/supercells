#include "slic_core.h"

/* Constructor. Nothing is done here. */
SlicCore::SlicCore() {

}

/* Destructor. Clear any present data. */
SlicCore::~SlicCore() {
  clear_data();
}

/* Clear the data as saved by the algorithm. */
void SlicCore::clear_data() {
  clusters.clear();
  distances.clear();
  centers.clear();
  centers_vals.clear();
  center_counts.clear();
  mat_dims.clear();
  new_clusters.clear();
  iter_mean_distance.clear();
  iter_diagnostics_enabled = false;
  vals_ptr = nullptr;
  bands = 0;
  dist_fn = nullptr;
  avg_fn = nullptr;
  avg_fun_name.clear();
}

/* Create centers of initial supercells, adds their values, etc. */
void SlicCore::create_centers(const std::vector<int>& mat_dims, const std::vector<double>& vals, int step) {
  for (int ncolcenter = step / 2; ncolcenter < mat_dims[1]; ncolcenter += step) {
    for (int nrowcenter = step / 2; nrowcenter < mat_dims[0]; nrowcenter += step) {
      std::vector<double> center; center.reserve(2);
      int ncell = ncolcenter + (nrowcenter * mat_dims[1]);
      std::vector<double> colour; colour.reserve(mat_dims[2]);
      for (int nval = 0; nval < mat_dims[2]; nval++) {
        double val = vals[ncell * mat_dims[2] + nval];
        colour.push_back(val);
      }

      std::vector<double> lm = find_local_minimum(vals, nrowcenter, ncolcenter);

      center.push_back(lm[0]);
      center.push_back(lm[1]);

      // To use the provided grid centers without snapping to local minima:
      // center[0] = nrowcenter;
      // center[1] = ncolcenter;

      centers.push_back(center);
      centers_vals.push_back(colour);
      center_counts.push_back(0);
    }
  }
}

/* Create centers of initial supercells based on the input provided by a user, adds their values, etc. */
void SlicCore::create_centers2(const std::vector<int>& mat_dims, const std::vector<double>& vals,
                              const std::vector<std::array<int, 2>>& input_centers) {
  for (size_t i = 0; i < input_centers.size(); i++) {
    int nrowcenter = input_centers[i][1];
    int ncolcenter = input_centers[i][0];
    std::vector<double> center; center.reserve(2);
    std::vector<double> colour; colour.reserve(mat_dims[2]);

    std::vector<double> lm = find_local_minimum(vals, nrowcenter, ncolcenter);

    int lm_row = static_cast<int>(lm[0]);
    int lm_col = static_cast<int>(lm[1]);
    int ncell = lm_col + (lm_row * mat_dims[1]);
    for (int nval = 0; nval < mat_dims[2]; nval++) {
      double val = vals[ncell * mat_dims[2] + nval];
      colour.push_back(val);
    }

    center.push_back(lm[0]);
    center.push_back(lm[1]);

    centers.push_back(center);
    centers_vals.push_back(colour);
    center_counts.push_back(0);
  }
}

/* Initialize the clusters, distances, and centers. */
void SlicCore::inits(const std::vector<int>& mat_dims_in, const std::vector<double>& vals,
                     const std::vector<std::array<int, 2>>& input_centers) {
  mat_dims = mat_dims_in;
  bands = mat_dims[2];
  vals_ptr = &vals;

  clusters.reserve(mat_dims[1]);
  distances.reserve(mat_dims[1]);
  for (int ncol = 0; ncol < mat_dims[1]; ncol++) {
    std::vector<int> cluster; cluster.reserve(mat_dims[0]);
    std::vector<double> distancemat; distancemat.reserve(mat_dims[0]);
    for (int nrow = 0; nrow < mat_dims[0]; nrow++) {
      cluster.push_back(-1);
      distancemat.push_back(FLT_MAX);
    }
    clusters.push_back(cluster);
    distances.push_back(distancemat);
  }

  if (input_centers.size() > 1) {
    create_centers2(mat_dims, vals, input_centers);
  } else {
    create_centers(mat_dims, vals, step);
  }
}

/* Compute the total (spatial and value) distance between a center and an individual pixel. */
double SlicCore::compute_dist(int ci, int y, int x, const std::vector<double>& values) const {
  double values_distance = dist_fn(centers_vals[ci], values);

  double y_dist = centers[ci][0] - y;
  double x_dist = centers[ci][1] - x;
  double spatial_distance = sqrt((y_dist * y_dist) + (x_dist * x_dist));

  double dist1 = values_distance / compactness;
  double dist2 = spatial_distance / step;

  return sqrt((dist1 * dist1) + (dist2 * dist2));
}

/* Find the pixel with the lowest gradient in a 3x3 surrounding. */
std::vector<double> SlicCore::find_local_minimum(const std::vector<double>& vals, int y, int x) {
  double min_grad = FLT_MAX;

  std::vector<double> loc_min(2);
  loc_min.at(0) = y;
  loc_min.at(1) = x;

  int rows = mat_dims[0];
  int cols = mat_dims[1];
  int bands = mat_dims[2];

  std::vector<double> colour1(bands);
  std::vector<double> colour2(bands);
  std::vector<double> colour3(bands);

  for (int i = x - 1; i < x + 2; i++) {
    for (int j = y - 1; j < y + 2; j++) {
      if (i < 0 || j < 0 || i + 1 >= cols || j + 1 >= rows) {
        continue;
      }

      int ncell1 = i + ((j + 1) * cols);
      int ncell2 = (i + 1) + (j * cols);
      int ncell3 = i + (j * cols);

      for (int nval = 0; nval < bands; nval++) {
        int idx1 = ncell1 * bands + nval;
        int idx2 = ncell2 * bands + nval;
        int idx3 = ncell3 * bands + nval;
        colour1[nval] = vals[idx1];
        colour2[nval] = vals[idx2];
        colour3[nval] = vals[idx3];
      }

      double new_grad = dist_fn(colour1, colour3) + dist_fn(colour2, colour3);

      if (new_grad < min_grad) {
        min_grad = new_grad;
        loc_min.at(0) = j;
        loc_min.at(1) = i;
      }
    }
  }
  return loc_min;
}
