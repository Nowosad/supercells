#include "slic_core.h"
#include "calc_stats.h"
#include <limits>

void SlicCore::generate_superpixels(const std::vector<int>& mat_dims_in, const double* vals_ptr_in,
                                    int step, double compactness, DistFn dist_fn_in, AvgFn avg_fn_in,
                                    const std::string& avg_fun_name_in, int iter,
                                    const std::vector<std::array<int, 2>>& input_centers,
                                    int verbose, bool diagnostics) {
  // Initialize state and seed cluster centers.
  if (verbose > 0) std::printf("Initialization: ");
  /* Clear previous data (if any), and re-initialize it. */
  clear_data();
  this->step = step;
  this->compactness = compactness;
  dist_fn = dist_fn_in;
  avg_fn = avg_fn_in;
  avg_fun_name = avg_fun_name_in;
  diagnostics_enabled = diagnostics;
  if (diagnostics_enabled) {
    iter_mean_distance.clear();
    iter_max_distance.clear();
    iter_frac_changed.clear();
  }
  inits(mat_dims_in, vals_ptr_in, input_centers);
  if (verbose > 0) std::printf("Completed\n");

  // Main SLIC loop: assign pixels -> update centers, for a fixed number of iterations.
  for (int itr = 0; itr < iter; itr++) {
    if (verbose > 0) std::printf("Iteration: %u/%u\r", itr + 1, iter);

    std::vector<std::vector<int> > prev_clusters;
    if (diagnostics_enabled) {
      prev_clusters = clusters;
    }

    // Reset per-pixel distances before assignment.
    for (int i = 0; i < mat_dims[1]; i++) {
      for (int j = 0; j < mat_dims[0]; j++) {
        distances[i][j] = FLT_MAX;
      }
    }

    // Assignment step: find the best center within a local window around each center.
    for (int l = 0; l < (int) centers.size(); l++) {
      /* Only compare to pixels in a 2 x step by 2 x step region. */
      for (int m = centers[l][1] - step; m < centers[l][1] + step; m++) {
        for (int n = centers[l][0] - step; n < centers[l][0] + step; n++) {
          if (m >= 0 && m < mat_dims[1] && n >= 0 && n < mat_dims[0]) {
            int ncell = m + (n * mat_dims[1]);

            std::vector<double> colour; colour.reserve(mat_dims[2]);

            int count_na = 0;
            for (int nval = 0; nval < mat_dims[2]; nval++) {
              double val = value_at(ncell, nval);
              colour.push_back(val);
              if (std::isnan(val)) {
                count_na += 1;
              }
            }
            /*check NAN*/
            if (count_na > 0) {
              continue;
            }
            double d = compute_dist(l, n, m, colour);

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

    // Update step: reset accumulators for center positions and values.
    for (int m = 0; m < (int) centers.size(); m++) {
      centers[m][0] = centers[m][1] = 0;
      for (int n = 0; n < (int) centers_vals[0].size(); n++) {
        centers_vals[m][n] = 0;
      }
      center_counts[m] = 0;
    }

    // Recompute centers using a non-mean (median/mean2/custom) aggregator.
    if (avg_fun_name != "mean") {
      IntToIntMap c_id_centers_vals;
      for (int l = 0; l < mat_dims[1]; l++) {
        for (int k = 0; k < mat_dims[0]; k++) {
          int c_id = clusters[l][k];
          if (c_id != -1) {
            int ncell = l + (k * mat_dims[1]);
            c_id_centers_vals.insert(std::make_pair(c_id, ncell));
            centers[c_id][0] += k;
            centers[c_id][1] += l;
            center_counts[c_id] += 1;
          }
        }
      }
      mapIter m_it, s_it;
      for (m_it = c_id_centers_vals.begin(); m_it != c_id_centers_vals.end(); m_it = s_it) {
        int c_id = (*m_it).first;
        std::pair<mapIter, mapIter> keyRange = c_id_centers_vals.equal_range(c_id);
        std::vector<std::vector<double> > centers_vals_c_id(mat_dims[2]);
        for (s_it = keyRange.first; s_it != keyRange.second; ++s_it) {
          int ncell = (*s_it).second;
          for (int nval = 0; nval < mat_dims[2]; nval++) {
            double val = value_at(ncell, nval);
            centers_vals_c_id[nval].push_back(val);
          }
        }
        for (int nval = 0; nval < mat_dims[2]; nval++) {
          // calculate
          if (avg_fun_name == "median") {
            centers_vals[c_id][nval] = median(centers_vals_c_id[nval]);
          } else if (avg_fun_name == "mean2") {
            centers_vals[c_id][nval] = mean(centers_vals_c_id[nval]);
          } else if (avg_fun_name.empty()) {
            centers_vals[c_id][nval] = avg_fn(centers_vals_c_id[nval]);
          }
        }
      }
      // /* Normalize the cluster centers. */
      for (int l = 0; l < (int) centers.size(); l++) {
        if (center_counts[l] > 0) {
          centers[l][0] /= center_counts[l];
          centers[l][1] /= center_counts[l];
        } else {
          centers[l][0] = -DBL_MAX;
          centers[l][1] = center_counts[l];
        }
      }
    } else if (avg_fun_name == "mean") {
      // Recompute centers using the standard mean aggregator.
      for (int l = 0; l < mat_dims[1]; l++) {
        for (int k = 0; k < mat_dims[0]; k++) {
          int c_id = clusters[l][k];

          if (c_id != -1) {
            int ncell = l + (k * mat_dims[1]);

            std::vector<double> colour; colour.reserve(mat_dims[2]);
            for (int nval = 0; nval < mat_dims[2]; nval++) {
              double val = value_at(ncell, nval);
              colour.push_back(val);
            }
            centers[c_id][0] += k;
            centers[c_id][1] += l;
            for (int nval = 0; nval < mat_dims[2]; nval++) {
              centers_vals[c_id][nval] += colour[nval];
            }
            center_counts[c_id] += 1;
          }
        }
      }
      /* Normalize the clusters (and their centers). */
      for (int l = 0; l < (int) centers.size(); l++) {
        if (center_counts[l] > 0) {
          centers[l][0] /= center_counts[l];
          centers[l][1] /= center_counts[l];
          for (int nval = 0; nval < mat_dims[2]; nval++) {
            centers_vals[l][nval] /= center_counts[l];
          }
        }
      }
    }

    // Diagnostics: track iteration-level convergence stats if enabled.
    if (diagnostics_enabled) {
      double sum_dist = 0.0;
      double max_dist = 0.0;
      int count = 0;
      int changed = 0;
      int total = 0;
      for (int i = 0; i < mat_dims[1]; i++) {
        for (int j = 0; j < mat_dims[0]; j++) {
          int curr = clusters[i][j];
          int prev = prev_clusters[i][j];
          if (curr != -1 || prev != -1) {
            total += 1;
            if (curr != prev) {
              changed += 1;
            }
          }
          if (curr != -1 && distances[i][j] != FLT_MAX) {
            sum_dist += distances[i][j];
            if (distances[i][j] > max_dist) {
              max_dist = distances[i][j];
            }
            count += 1;
          }
        }
      }
      if (count > 0) {
        iter_mean_distance.push_back(sum_dist / count);
        iter_max_distance.push_back(max_dist);
      } else {
        double na = std::numeric_limits<double>::quiet_NaN();
        iter_mean_distance.push_back(na);
        iter_max_distance.push_back(na);
      }
      if (total > 0) {
        iter_frac_changed.push_back((double) changed / total);
      } else {
        iter_frac_changed.push_back(std::numeric_limits<double>::quiet_NaN());
      }
    }
  }
  if (verbose > 0) std::printf("\n");
}
