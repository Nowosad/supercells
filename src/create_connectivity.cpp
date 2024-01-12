#include "slic.h"

void Slic::create_connectivity(doubles_matrix<> vals, cpp11::function avg_fun_fun, std::string& avg_fun_name, int minarea, int verbose) {
  if (verbose > 0) Rprintf("Cleaning connectivity: ");
  int label = 0;
  int adjlabel = -1;
  const int dx4[4] = {-1,  0,  1,  0};
  const int dy4[4] = { 0, -1,  0,  1};

  if (minarea == 0) {
    minarea = (mat_dims[1] * mat_dims[0]) / ((int)centers.size());
    minarea = minarea >> 2;
  }

  for (int i = 0; i < mat_dims[1]; i++) {
    vector<int> ncl; ncl.reserve(mat_dims[0]);
    for (int j = 0; j < mat_dims[0]; j++) {
      ncl.push_back(-1);
    }
    new_clusters.push_back(ncl);
  }

  for (int i = 0; i < mat_dims[1]; i++) {
    for (int j = 0; j < mat_dims[0]; j++) {

      if (new_clusters[i][j] == -1 && clusters[i][j] != -1) {

        new_clusters[i][j] = label;

        vector<int> element(2);
        element.at(0) = j;
        element.at(1) = i;

        vector<vector<int> > elements;
        elements.push_back(element);

        /* Find an adjacent label, for possible use later. */
        for (int k = 0; k < 4; k++) {
          int x = elements[0][1] + dx4[k], y = elements[0][0] + dy4[k];

          if (x >= 0 && x < mat_dims[1] && y >= 0 && y < mat_dims[0]) {
            if (new_clusters[x][y] >= 0) {
              adjlabel = new_clusters[x][y];
            // } else {
            //   adjlabel = -1;
            }
          }
        }

        int count = 1;
        for (int c = 0; c < count; c++) {
          for (int k = 0; k < 4; k++) {
            int x = elements[c][1] + dx4[k], y = elements[c][0] + dy4[k];

            if (x >= 0 && x < mat_dims[1] && y >= 0 && y < mat_dims[0]) {
              if (new_clusters[x][y] == -1 && clusters[i][j] == clusters[x][y]) {
                vector<int> element2(2);
                element2.at(0) = y;
                element2.at(1) = x;

                elements.push_back(element2);
                new_clusters[x][y] = label;
                count += 1;
              }
            }
          }
        }

        /* Use the earlier found adjacent label if a segment size is
         smaller than a limit. */
        if (count <= minarea) {
          for (int c = 0; c < count; c++) {
            new_clusters[elements[c][1]][elements[c][0]] = adjlabel;
          }
          label = label - 1;
        }
        label = label + 1;
        adjlabel = -1;
      }
    }
  }

  clusters = new_clusters;

  vector<vector<int> > new_centers;
  vector<vector<double> > new_centers_vals;
  vector<int> new_center_counts(label);

  /* Clear the center values. */
  /* Clear the center _vals values. */
  for (int m = 0; m < (int) label; m++) {

    vector<int> new_center(2);
    new_center[0] = new_center[1] = 0;
    new_centers.push_back(new_center);

    vector<double> new_center_val;
    for (int n = 0; n < (int) mat_dims[2]; n++){
      new_center_val.push_back(0);
    }
    new_centers_vals.push_back(new_center_val);
    new_center_counts[m] = 0;
  }

  if (avg_fun_name != "mean"){
    multimap <int, int> new_c_id_centers_vals;
    for (int l = 0; l < (int) new_clusters.size(); l++) {
      for (int k = 0; k < (int) new_clusters[0].size(); k++) {
        int c_id = new_clusters[l][k];
        if (c_id != -1) {
          int ncell = l + (k * mat_dims[1]);
          new_c_id_centers_vals.insert(make_pair(c_id, ncell));
          new_centers[c_id][0] += k;
          new_centers[c_id][1] += l;
          new_center_counts[c_id] += 1;
        }
      }
    }
    mapIter m_it, s_it;
    for (m_it = new_c_id_centers_vals.begin();  m_it != new_c_id_centers_vals.end();  m_it = s_it){
      int c_id = (*m_it).first;
      pair<mapIter, mapIter> keyRange = new_c_id_centers_vals.equal_range(c_id);
      vector<vector<double> > new_c_id_centers_vals(mat_dims[2]);
      // Iterate over all map elements with key == theKey
      for (s_it = keyRange.first;  s_it != keyRange.second;  ++s_it){
        int ncell = (*s_it).second;
        for (int nval = 0; nval < mat_dims[2]; nval++){
          double val = vals(ncell, nval);
          new_c_id_centers_vals[nval].push_back(val);
        }
      }
      for (int nval = 0; nval < mat_dims[2]; nval++){
        // calculate
        if (avg_fun_name == "median"){
          new_centers_vals[c_id][nval] = median(new_c_id_centers_vals[nval]);
        } else if (avg_fun_name == "mean2"){
          new_centers_vals[c_id][nval] = mean(new_c_id_centers_vals[nval]);
        } else if (avg_fun_name.empty()){
          // use user-defined function
          new_centers_vals[c_id][nval] = avg_fun_fun(new_c_id_centers_vals[nval]);
        }
      }
    }

    // /* Normalize the clusters. */
    for (int l = 0; l < (int) label; l++) {
      new_centers[l][0] /= new_center_counts[l];
      new_centers[l][1] /= new_center_counts[l];
    }
  } else if (avg_fun_name == "mean"){
    // /* Compute the new cluster centers. */
    for (int l = 0; l < (int) new_clusters.size(); l++) {
      for (int k = 0; k < (int) new_clusters[0].size(); k++) {
        int c_id = new_clusters[l][k];

        if (c_id != -1) {
          int ncell = l + (k * mat_dims[1]);

          vector<double> colour;
          for (int nval = 0; nval < mat_dims[2]; nval++){
            double val = vals(ncell, nval);
            colour.push_back(val);
          }

          new_centers[c_id][0] += k;
          new_centers[c_id][1] += l;
          for (int nval = 0; nval < mat_dims[2]; nval++){
            new_centers_vals[c_id][nval] += colour[nval];
          }

          new_center_counts[c_id] += 1;
        }
      }
    }
    // /* Normalize the clusters. */
    for (int l = 0; l < (int) label; l++) {
      new_centers[l][0] /= new_center_counts[l];
      new_centers[l][1] /= new_center_counts[l];
      for (int nval = 0; nval < mat_dims[2]; nval++){
        new_centers_vals[l][nval] /= new_center_counts[l];
      }
    }
  }

  centers = new_centers;
  centers_vals = new_centers_vals;
  if (verbose > 0) Rprintf("Completed\n");
}
