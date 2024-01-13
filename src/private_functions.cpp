#include "slic.h"

/* Constructor. Nothing is done here. */
Slic::Slic() {

}

/* Destructor. Clear any present data. */
Slic::~Slic() {
  clear_data();
}

/* Clear the data as saved by the algorithm. */
void Slic::clear_data() {
  clusters.clear();
  distances.clear();
  centers.clear();
  centers_vals.clear();
  center_counts.clear();
  mat_dims.clear();
}

void Slic::create_centers(vector<int> mat_dims, doubles_matrix<> vals,
                          std::string& dist_name, cpp11::function dist_fun, double step) {
  for (int ncolcenter = step/2; ncolcenter < mat_dims[1]; ncolcenter += step){
    for (int nrowcenter = step/2; nrowcenter < mat_dims[0]; nrowcenter += step){
      vector<int> center; center.reserve(2);
      int ncell = ncolcenter + (nrowcenter * mat_dims[1]);
      vector<double> colour; colour.reserve(mat_dims[2]);
      for (int nval = 0; nval < mat_dims[2]; nval++){
        double val = vals(ncell, nval);
        colour.push_back(val);
      }

      vector<int> lm = find_local_minimum(vals, nrowcenter, ncolcenter, dist_name, dist_fun);

      /* Generate the center vector. */
      center.push_back(lm[0]);
      center.push_back(lm[1]);

      /* Append to vector of centers. */
      centers.push_back(center);
      centers_vals.push_back(colour);
      center_counts.push_back(0);
    }
  }
}

void Slic::create_centers2(vector<int> mat_dims,
                           doubles_matrix<> vals, std::string& dist_name,
                           cpp11::function dist_fun,
                           integers_matrix<> input_centers) {
  int ncell = 0;
  for (int i = 0; i < input_centers.nrow(); i++){
    int nrowcenter = input_centers(i, 1); int ncolcenter = input_centers(i, 0);
    vector<int> center; center.reserve(2);
    vector<double> colour; colour.reserve(mat_dims[2]);
    for (int nval = 0; nval < mat_dims[2]; nval++){
      double val = vals(ncell, nval);
      colour.push_back(val);
    }
    ncell++;

    vector<int> lm = find_local_minimum(vals, nrowcenter, ncolcenter, dist_name, dist_fun);

    /* Generate the center vector. */
    center.push_back(lm[0]);
    center.push_back(lm[1]);

    centers.push_back(center);
    centers_vals.push_back(colour);
    center_counts.push_back(0);
  }
}

void Slic::inits(integers mat, doubles_matrix<> vals,
                 std::string& dist_name, cpp11::function dist_fun,
                 integers_matrix<> input_centers) {
  // cout << "inits" << endl;

  mat_dims.reserve(3);
  mat_dims.push_back(mat.at(0));
  mat_dims.push_back(mat.at(1));
  mat_dims.push_back(vals.ncol());

  clusters.reserve(mat_dims[1]);
  distances.reserve(mat_dims[1]);
  for (int ncol = 0; ncol < mat_dims[1]; ncol++){
    vector<int> cluster; cluster.reserve(mat_dims[0]);
    vector<double> distancemat; distancemat.reserve(mat_dims[0]);
    for (int nrow = 0; nrow < mat_dims[0]; nrow++){
      cluster.push_back(-1);
      distancemat.push_back(FLT_MAX);
    }
    clusters.push_back(cluster);
    distances.push_back(distancemat);
  }

  if (input_centers.nrow() > 1){
    Slic::create_centers2(mat_dims, vals, dist_name, dist_fun, input_centers);
  } else{
    Slic::create_centers(mat_dims, vals, dist_name, dist_fun, step);
  }
}

double Slic::compute_dist(int& ci, int& y, int& x, vector<double>& values,
                          std::string& dist_name, cpp11::function dist_fun) {

  /*vals distance*/
  double values_distance = get_vals_dist(centers_vals[ci], values, dist_name, dist_fun);

  /*coords distance*/
  int y_dist = centers[ci][0] - y;
  int x_dist = centers[ci][1] - x;
  double spatial_distance = sqrt((y_dist * y_dist) + (x_dist * x_dist));

  double dist1 = values_distance / compactness;
  double dist2 = spatial_distance / step;

  return sqrt((dist1 * dist1) + (dist2 * dist2));
}

vector<int> Slic::find_local_minimum(doubles_matrix<> vals, int& y, int& x,
                                     std::string& dist_name, cpp11::function dist_fun) {
  double min_grad = FLT_MAX;

  vector<int> loc_min(2);
  loc_min.at(0) = y;
  loc_min.at(1) = x;

  for (int i = x - 1; i < x + 2; i++) {
    for (int j = y - 1; j < y + 2; j++) {

      int ncell1 = i + ((j + 1) * mat_dims[1]);
      int ncell2 = (i + 1) + (j * mat_dims[1]);
      int ncell3 = i + (j * mat_dims[1]);

      vector<double> colour1; colour1.reserve(mat_dims[2]);
      vector<double> colour2; colour2.reserve(mat_dims[2]);
      vector<double> colour3; colour3.reserve(mat_dims[2]);

      if (ncell1 < vals.nrow() && ncell2 < vals.nrow() && ncell3 < vals.nrow()){
        for (int nval = 0; nval < mat_dims[2]; nval++){
          double val1 = vals(ncell1, nval);
          double val2 = vals(ncell2, nval);
          double val3 = vals(ncell3, nval);
          colour1.push_back(val1);
          colour2.push_back(val2);
          colour3.push_back(val3);
        }

        /* Compute horizontal and vertical gradients and keep track of the
         minimum. */
        double new_grad = get_vals_dist(colour1, colour3, dist_name, dist_fun) + get_vals_dist(colour2, colour3, dist_name, dist_fun);

        if (new_grad < min_grad) {
          min_grad = new_grad;
          loc_min.at(0) = j;
          loc_min.at(1) = i;
        }
      }
    }
  }
  return loc_min;
}
