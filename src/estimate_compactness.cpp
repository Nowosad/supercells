#include "slic.h"

writable::doubles Slic::estimate_compactness(integers mat, doubles_matrix<> vals, double step,
                                std::string& dist_name, cpp11::function dist_fun,
                                integers_matrix<> input_centers){

  this->step = step;

  /* Clear previous data (if any), and re-initialize it. */
  clear_data();
  inits(mat, vals, dist_name, dist_fun, input_centers);

  /* Reset distance values. */
  for (int i = 0; i < mat_dims[1]; i++) {
    for (int j = 0; j < mat_dims[0]; j++) {
      distances[i][j] = FLT_MAX;
    }
  }
  /////////////////////////////////
  //create vector that stores maximum value distance for each center
  vector<double> max_dist; max_dist.reserve(centers.size());

  for (int l = 0; l < (int) centers.size(); l++) {
    double max_d = 0; //initialize maximum distance as 0

    /* Only compare to pixels in a 2 x step by 2 x step region. */
    for (int m = centers[l][1] - step; m < centers[l][1] + step; m++) {
      for (int n = centers[l][0] - step; n < centers[l][0] + step; n++) {

        if (m >= 0 && m < mat_dims[1] && n >= 0 && n < mat_dims[0]) {

          int ncell = m + (n * mat_dims[1]);

          vector<double> colour; colour.reserve(mat_dims[2]);

          int count_na = 0;
          for (int nval = 0; nval < mat_dims[2]; nval++){
            double val = vals(ncell, nval);
            colour.push_back(val);
            int nanr = is_na(val);
            count_na = count_na + nanr;
          }
          /*check NAN*/
          if (count_na > 0){
            continue;
          }
          double d = get_vals_dist(centers_vals[l], colour, dist_name, dist_fun);
          if(d > max_d){
            max_d = d; //update maximum distance for a given center if value is greater than previous maximum
          }
        }
      }
      max_dist[l] = max_d; //store maximum value for corresponding center
    }
  }
  //summarize across cluster centers
  writable::doubles est_compactness(4);
  est_compactness[0] = *std::min_element(max_dist.begin(), max_dist.end()); //minimum across cluster centers
  est_compactness[1] = median(max_dist); //median across cluster centers
  est_compactness[2] = mean(max_dist); //mean across cluster centers
  est_compactness[3] = *std::max_element(max_dist.begin(), max_dist.end()); //maximum across cluster centers
  writable::strings names(4); //set names of output vector
  names[0] = "min";
  names[1] = "median";
  names[2] = "mean";
  names[3] = "max";
  est_compactness.attr("names") = names;
  return est_compactness; //return estimated compactness summary
}







