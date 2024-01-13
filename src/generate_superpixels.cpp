#include "slic.h"

void Slic::generate_superpixels(integers mat, doubles_matrix<> vals, double step, double compactness,
                                std::string& dist_name, cpp11::function dist_fun,
                                cpp11::function avg_fun_fun, std::string& avg_fun_name, int iter,
                                integers_matrix<> input_centers,
                                int verbose){
  this->step = step;
  this->compactness = compactness;

  if (verbose > 0) Rprintf("Initialization: ");
  /* Clear previous data (if any), and re-initialize it. */
  clear_data();
  inits(mat, vals, dist_name, dist_fun, input_centers);
  if (verbose > 0) Rprintf("Completed\n");

  /* Run for provided number of iterations. */
  for (int itr = 0; itr < iter; itr++) {

    if (verbose > 0) Rprintf("Iteration: %u/%u\r", itr + 1, iter);

    /* Reset distance values. */
    for (int i = 0; i < mat_dims[1]; i++) {
      for (int j = 0; j < mat_dims[0]; j++) {
        distances[i][j] = FLT_MAX;
      }
    }

    for (int l = 0; l < (int) centers.size(); l++) {
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
            double d = compute_dist(l, n, m, colour, dist_name, dist_fun);

            /* Update cluster allocation if the cluster minimizes the
             distance. */
            if (d < distances[m][n]) {
              distances[m][n] = d;
              clusters[m][n] = l;
            }
          }
        }
      }
    }

    /* Clear the center values. */
    /* Clear the center_vals values. */
    for (int m = 0; m < (int) centers.size(); m++) {
      centers[m][0] = centers[m][1] = 0;
      for (int n = 0; n < (int) centers_vals[0].size(); n++){
        centers_vals[m][n] = 0;
      }
      center_counts[m] = 0;
    }
    
    /* Compute the new cluster centers when the average function is not mean. */
    if (avg_fun_name != "mean"){
      multimap <int, int> c_id_centers_vals;
      for (int l = 0; l < mat_dims[1]; l++) {
        for (int k = 0; k < mat_dims[0]; k++) {
          int c_id = clusters[l][k];
          if (c_id != -1) {
            int ncell = l + (k * mat_dims[1]);
            c_id_centers_vals.insert(make_pair(c_id, ncell));
            centers[c_id][0] += k;
            centers[c_id][1] += l;
            center_counts[c_id] += 1;
          }
        }
      }
      mapIter m_it, s_it;
      for (m_it = c_id_centers_vals.begin();  m_it != c_id_centers_vals.end();  m_it = s_it){
        int c_id = (*m_it).first;
        pair<mapIter, mapIter> keyRange = c_id_centers_vals.equal_range(c_id);
        vector<vector<double> > centers_vals_c_id(mat_dims[2]);
        // Iterate over all map elements with key == theKey
        for (s_it = keyRange.first; s_it != keyRange.second; ++s_it){
          int ncell = (*s_it).second;
            for (int nval = 0; nval < mat_dims[2]; nval++){
              double val = vals(ncell, nval);
              centers_vals_c_id[nval].push_back(val);
            }
        }
        for (int nval = 0; nval < mat_dims[2]; nval++){
          // calculate
          if (avg_fun_name == "median"){
            centers_vals[c_id][nval] = median(centers_vals_c_id[nval]);
          } else if (avg_fun_name == "mean2"){
            centers_vals[c_id][nval] = mean(centers_vals_c_id[nval]);
          } else if (avg_fun_name.empty()){
            centers_vals[c_id][nval] = avg_fun_fun(centers_vals_c_id[nval]);
          }
        }
      }
      // /* Normalize the cluster centers. */
      for (int l = 0; l < (int) centers.size(); l++) {
        if (center_counts[l] > 0){
          centers[l][0] /= center_counts[l];
          centers[l][1] /= center_counts[l];
        } else {
          centers[l][0] = INT_MIN;
          centers[l][1] = center_counts[l];
        }
      }
    } else if (avg_fun_name == "mean"){
      /* Compute the new cluster centers when the average function is mean. */
      for (int l = 0; l < mat_dims[1]; l++) {
        for (int k = 0; k < mat_dims[0]; k++) {
          int c_id = clusters[l][k];

          if (c_id != -1) {
            int ncell = l + (k * mat_dims[1]);

            vector<double> colour; colour.reserve(mat_dims[2]);
            for (int nval = 0; nval < mat_dims[2]; nval++){
              double val = vals(ncell, nval);
              colour.push_back(val);
            }
            centers[c_id][0] += k;
            centers[c_id][1] += l;
            for (int nval = 0; nval < mat_dims[2]; nval++){
              centers_vals[c_id][nval] += colour[nval];
            }
            center_counts[c_id] += 1;
          }
        }
      }
      /* Normalize the clusters (and their centers). */
      for (int l = 0; l < (int) centers.size(); l++) {
        if (center_counts[l] > 0){
          centers[l][0] /= center_counts[l];
          centers[l][1] /= center_counts[l];
          for (int nval = 0; nval < mat_dims[2]; nval++){
            centers_vals[l][nval] /= center_counts[l];
          }
        }
      }
    }
  }
  if (verbose > 0) Rprintf("\n");
}
