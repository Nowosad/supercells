#include "slic_core.h"
#include "calc_stats.h"

void SlicCore::create_connectivity(const std::vector<double>& vals, AvgFn avg_fn_in,
                                   const std::string& avg_fun_name_in, int minarea, int verbose) {
  avg_fn = avg_fn_in;
  avg_fun_name = avg_fun_name_in;
  if (verbose > 0) std::printf("Cleaning connectivity: ");
  int label = 0;
  int adjlabel = -1;
  const int dx4[4] = {-1,  0,  1,  0};
  const int dy4[4] = { 0, -1,  0,  1};

  if (centers.empty()) {
    if (verbose > 0) std::printf("No centers available\n");
    return;
  }

  if (minarea == 0) {
    minarea = (mat_dims[1] * mat_dims[0]) / ((int)centers.size());
    minarea = minarea >> 2;
  }

  new_clusters.assign(mat_dims[1], std::vector<int>(mat_dims[0], -1));

  std::vector<std::array<int, 2>> elements;
  for (int i = 0; i < mat_dims[1]; i++) {
    for (int j = 0; j < mat_dims[0]; j++) {

      if (new_clusters[i][j] == -1 && clusters[i][j] != -1) {

        new_clusters[i][j] = label;

        elements.clear();
        elements.push_back({j, i});

        /* Find an adjacent label, for possible use later. */
        for (int k = 0; k < 4; k++) {
          int x = elements[0][1] + dx4[k], y = elements[0][0] + dy4[k];

          if (x >= 0 && x < mat_dims[1] && y >= 0 && y < mat_dims[0]) {
            if (new_clusters[x][y] >= 0) {
              adjlabel = new_clusters[x][y];
            }
          }
        }

        int count = 1;
        for (int c = 0; c < count; c++) {
          for (int k = 0; k < 4; k++) {
            int x = elements[c][1] + dx4[k], y = elements[c][0] + dy4[k];

            if (x >= 0 && x < mat_dims[1] && y >= 0 && y < mat_dims[0]) {
              if (new_clusters[x][y] == -1 && clusters[i][j] == clusters[x][y]) {
                elements.push_back({y, x});
                new_clusters[x][y] = label;
                count += 1;
              }
            }
          }
        }

        /* Use the earlier found adjacent label if a segment size is
         smaller than a limit. */
        if (count <= minarea) {
          int target = adjlabel;
          if (target == -1) {
            target = label;
          }
          for (int c = 0; c < count; c++) {
            new_clusters[elements[c][1]][elements[c][0]] = target;
          }
          if (adjlabel != -1) {
            label = label - 1;
          }
        }
        label = label + 1;
        adjlabel = -1;
      }
    }
  }

  clusters = new_clusters;

  const auto prev_centers = centers;
  const auto prev_centers_vals = centers_vals;
  std::vector<std::vector<double>> new_centers(label, std::vector<double>(2, 0));
  std::vector<std::vector<double>> new_centers_vals(label, std::vector<double>(mat_dims[2], 0));
  std::vector<int> new_center_counts(label, 0);

  std::vector<std::vector<int>> cluster_cells(label);
  std::vector<std::vector<double>> new_c_id_centers_vals_buf(mat_dims[2]);
  std::vector<double> colour(mat_dims[2]);
  if (avg_fun_name != "mean") {
    for (int l = 0; l < (int) new_clusters.size(); l++) {
      for (int k = 0; k < (int) new_clusters[0].size(); k++) {
        int c_id = new_clusters[l][k];
        if (c_id != -1) {
          int ncell = l + (k * mat_dims[1]);
          cluster_cells[c_id].push_back(ncell);
          new_centers[c_id][0] += k;
          new_centers[c_id][1] += l;
          new_center_counts[c_id] += 1;
        }
      }
    }
    for (int c_id = 0; c_id < label; c_id++) {
      if (cluster_cells[c_id].empty()) {
        continue;
      }
      for (int nval = 0; nval < mat_dims[2]; nval++) {
        new_c_id_centers_vals_buf[nval].clear();
      }
      for (int ncell : cluster_cells[c_id]) {
        for (int nval = 0; nval < mat_dims[2]; nval++) {
            double val = vals[ncell * mat_dims[2] + nval];
            new_c_id_centers_vals_buf[nval].push_back(val);
        }
      }
      for (int nval = 0; nval < mat_dims[2]; nval++) {
        // calculate
        if (avg_fun_name == "median") {
          new_centers_vals[c_id][nval] = median(new_c_id_centers_vals_buf[nval]);
        } else if (avg_fun_name == "mean2") {
          new_centers_vals[c_id][nval] = mean(new_c_id_centers_vals_buf[nval]);
        } else if (avg_fun_name.empty()) {
          new_centers_vals[c_id][nval] = avg_fn(new_c_id_centers_vals_buf[nval]);
        }
      }
    }

    // /* Normalize the clusters. */
    for (int l = 0; l < (int) label; l++) {
      if (new_center_counts[l] > 0) {
        new_centers[l][0] /= new_center_counts[l];
        new_centers[l][1] /= new_center_counts[l];
      } else if (l < (int) prev_centers.size()) {
        new_centers[l] = prev_centers[l];
        new_centers_vals[l] = prev_centers_vals[l];
      }
    }
  } else if (avg_fun_name == "mean") {
    // /* Compute the new cluster centers. */
    for (int l = 0; l < (int) new_clusters.size(); l++) {
      for (int k = 0; k < (int) new_clusters[0].size(); k++) {
        int c_id = new_clusters[l][k];

        if (c_id != -1) {
          int ncell = l + (k * mat_dims[1]);

          for (int nval = 0; nval < mat_dims[2]; nval++) {
            double val = vals[ncell * mat_dims[2] + nval];
            colour[nval] = val;
          }

          new_centers[c_id][0] += k;
          new_centers[c_id][1] += l;
          for (int nval = 0; nval < mat_dims[2]; nval++) {
            new_centers_vals[c_id][nval] += colour[nval];
          }

          new_center_counts[c_id] += 1;
        }
      }
    }
    // /* Normalize the clusters. */
    for (int l = 0; l < (int) label; l++) {
      if (new_center_counts[l] > 0) {
        new_centers[l][0] /= new_center_counts[l];
        new_centers[l][1] /= new_center_counts[l];
        for (int nval = 0; nval < mat_dims[2]; nval++) {
          new_centers_vals[l][nval] /= new_center_counts[l];
        }
      } else if (l < (int) prev_centers.size()) {
        new_centers[l] = prev_centers[l];
        new_centers_vals[l] = prev_centers_vals[l];
      }
    }
  }

  centers = new_centers;
  centers_vals = new_centers_vals;
  if (verbose > 0) std::printf("Completed\n");
}
