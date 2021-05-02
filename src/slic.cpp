#include "slic.h"

/*
 * Constructor. Nothing is done here.
 */
Slic::Slic() {

}

/*
 * Destructor. Clear any present data.
 */
Slic::~Slic() {
  clear_data();
}

/*
 * Clear the data as saved by the algorithm.
 *
 * Input : -
 * Output: -
 */
void Slic::clear_data() {
  clusters.clear();
  distances.clear();
  centers.clear();
  centers_vals.clear();
  center_counts.clear();
}

void Slic::inits(integers_matrix m, integers_matrix vals){
  // cout << "inits" << endl;

  for (int ncol = 0; ncol < m.ncol(); ncol++){
    vector<int> cluster;
    vector<double> distancemat;
    for (int nrow = 0; nrow < m.nrow(); nrow++){
      cluster.push_back(-1);
      distancemat.push_back(FLT_MAX);
    }
    clusters.push_back(cluster);
    distances.push_back(distancemat);
  }

  for (int ncolcenter = step; ncolcenter < m.ncol() - step/2; ncolcenter += step){
    for (int nrowcenter = step; nrowcenter < m.nrow() - step/2; nrowcenter += step){
      vector<double> center;
      // cout << ncolcenter << " " << nrowcenter << endl;

      int ncell = ncolcenter + (nrowcenter * m.ncol());

      vector<double> colour;
      for (int nval = 0; nval < vals.ncol(); nval++){
        double val = vals(ncell, nval);
        colour.push_back(val);
      }

      vector<int> lm = find_local_minimum(vals, m, nrowcenter, ncolcenter);

      /* Generate the center vector. */
      center.push_back(lm[0]);
      center.push_back(lm[1]);
      // center.push_back(colour);

      /* Generate the center vector. */
      // center.push_back(nrowcenter);
      // center.push_back(ncolcenter);
      // center.insert(center.end(), colour.begin(), colour.end());

      /* Append to vector of centers. */
      centers.push_back(center);
      centers_vals.push_back(colour);
      center_counts.push_back(0);
    }
  }
}

double Slic::euclidean(vector<double> values1, vector<double> values2){

    int len1 = values1.size();
    // int len2 = values2.size();
    double dist = 0.0;
    double diff = 0.0;

    for (int i = 0; i < len1; i++){
        diff = values1[i] - values2[i];
        dist += pow(diff, 2);
    }
    return sqrt(dist);
}

double Slic::compute_dist(int ci, int y, int x, vector<double> values) {

  double dc = euclidean(centers_vals[ci], values);
  double ds = sqrt(pow(centers[ci][0] - y, 2) + pow(centers[ci][1] - x, 2));
  // cout << "distance: " << dc << endl;

  return sqrt(pow(dc / nc, 2) + pow(ds / ns, 2));
}

vector<int> Slic::find_local_minimum(integers_matrix vals, integers_matrix mat, int y, int x) {
  double min_grad = FLT_MAX;
  // int min_grad = -1;

  vector<int> loc_min(2);
  loc_min.at(0) = y;
  loc_min.at(1) = x;

  // cout << "loc min start: " << loc_min[0] << " " << loc_min[1] << endl;

  for (int i = x - 1; i < x + 2; i++) {
    for (int j = y - 1; j < y + 2; j++) {

      int ncell1 = i + ((j + 1) * mat.ncol());
      int ncell2 = (i + 1) + (j * mat.ncol());
      int ncell3 = i + (j * mat.ncol());

      vector<double> colour1;
      vector<double> colour2;
      vector<double> colour3;
      for (int nval = 0; nval < vals.ncol(); nval++){
        double val1 = vals(ncell1, nval);
        double val2 = vals(ncell2, nval);
        double val3 = vals(ncell3, nval);
        colour1.push_back(val1);
        colour2.push_back(val2);
        colour3.push_back(val3);
      }

      /* Compute horizontal and vertical gradients and keep track of the
       minimum. */
      double new_grad = euclidean(colour1, colour3) + euclidean(colour2, colour3);
      // cout << "new_grad: " << new_grad << " min_grad: " << min_grad << endl;

      if (new_grad < min_grad) {
        min_grad = new_grad;
        loc_min.at(0) = j;
        loc_min.at(1) = i;
      }
    }
  }
  // cout << "loc min end: " << loc_min[0] << " " << loc_min[1] << endl;

  return loc_min;
}

void Slic::generate_superpixels(integers_matrix mat, integers_matrix vals, double step, int nc){
  // cout << "generate_superpixels" << endl;
  this->step = step;
  this->nc = nc;
  this->ns = step;

  /* Clear previous data (if any), and re-initialize it. */
  clear_data();
  inits(mat, vals);

  /* Run EM for 10 iterations (as prescribed by the algorithm). */
  for (int iter = 0; iter < NR_ITERATIONS; iter++) {
    /* Reset distance values. */
    for (int i = 0; i < mat.ncol(); i++) {
      for (int j = 0; j < mat.nrow(); j++) {
        distances[i][j] = FLT_MAX;
      }
    }
    for (int l = 0; l < (int) centers.size(); l++) {
      /* Only compare to pixels in a 2 x step by 2 x step region. */
      for (int m = centers[l][1] - step; m < centers[l][1] + step; m++) {
        for (int n = centers[l][0] - step; n < centers[l][0] + step; n++) {

          if (m >= 0 && m < mat.ncol() && n >= 0 && n < mat.nrow()) {
            // int colour = mat(n, m);

            int ncell = m + (n * mat.ncol());

            vector<double> colour;
            for (int nval = 0; nval < vals.ncol(); nval++){
              double val = vals(ncell, nval);
              colour.push_back(val);
            }

            double d = compute_dist(l, n, m, colour);

            // cout << "distance: " << d << endl;

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
    /* Clear the center _vals values. */
    for (int m = 0; m < (int) centers.size(); m++) {
      centers[m][0] = centers[m][1] = 0;
      for (int n = 0; n < (int) centers_vals[0].size(); n++){
        centers_vals[m][n] = 0;
      }
      center_counts[m] = 0;
    }

    // /* Compute the new cluster centers. */
    for (int l = 0; l < mat.ncol(); l++) {
      for (int k = 0; k < mat.nrow(); k++) {
        int c_id = clusters[l][k];

        if (c_id != -1) {
          int ncell = l + (k * mat.ncol());

          vector<double> colour;
          for (int nval = 0; nval < vals.ncol(); nval++){
            double val = vals(ncell, nval);
            colour.push_back(val);
          }
          centers[c_id][0] += k;
          centers[c_id][1] += l;
          for (int nval = 0; nval < vals.ncol(); nval++){
            centers_vals[c_id][nval] += colour[nval];
          }
          center_counts[c_id] += 1;
        }
      }
    }
    // /* Normalize the clusters. */
    for (int l = 0; l < (int) centers.size(); l++) {
      centers[l][0] /= center_counts[l];
      centers[l][1] /= center_counts[l];
      for (int nval = 0; nval < vals.ncol(); nval++){
        centers_vals[l][nval] /= center_counts[l];
      }
    }
  }
}

void Slic::create_connectivity(integers_matrix mat, integers_matrix vals) {
  int label = 0;
  int adjlabel = 0;
  const int lims = (mat.ncol() * mat.nrow()) / ((int)centers.size());
  const int dx4[4] = {-1,  0,  1,  0};
  const int dy4[4] = { 0, -1,  0,  1};

  for (int i = 0; i < mat.ncol(); i++) {
    vector<int> nc;
    for (int j = 0; j < mat.nrow(); j++) {
      nc.push_back(-1);
    }
    new_clusters.push_back(nc);
  }

  for (int i = 0; i < mat.ncol(); i++) {
    for (int j = 0; j < mat.nrow(); j++) {

      if (new_clusters[i][j] == -1) {

        new_clusters[i][j] = label;

        vector<vector<int> > elements;
        vector<int> element;
        element.push_back(j);
        element.push_back(i);
        elements.push_back(element);

        /* Find an adjacent label, for possible use later. */
        for (int k = 0; k < 4; k++) {
          int x = elements[0][1] + dx4[k], y = elements[0][0] + dy4[k];

          if (x >= 0 && x < mat.ncol() && y >= 0 && y < mat.nrow()) {
            if (new_clusters[x][y] >= 0) {
              adjlabel = new_clusters[x][y];
              // cout << "adjlabel" << adjlabel <<endl;

            }
          }
        }

        int count = 1;
        for (int c = 0; c < count; c++) {
          for (int k = 0; k < 4; k++) {
            int x = elements[c][1] + dx4[k], y = elements[c][0] + dy4[k];

            if (x >= 0 && x < mat.ncol() && y >= 0 && y < mat.nrow()) {
              if (new_clusters[x][y] == -1 && clusters[i][j] == clusters[x][y]) {
                vector<int> element2;
                element2.push_back(y);
                element2.push_back(x);
                elements.push_back(element2);
                new_clusters[x][y] = label;
                count += 1;
              }
            }
          }
        }

        /* Use the earlier found adjacent label if a segment size is
         smaller than a limit. */
        if (count <= lims >> 2) {
          for (int c = 0; c < count; c++) {
            new_clusters[elements[c][1]][elements[c][0]] = adjlabel;
          }
          label = label - 1;
        }
        label = label + 1;
        // cout << "label" << label << "";

      }
    }
  }

  // cout << "clusters" << clusters.size() << " " << clusters[0].size() <<endl;
  // cout << "new_clusters" << new_clusters.size() << " " << new_clusters[0].size() <<endl;

  clusters = new_clusters;

  /* Clear the center values. */
  /* Clear the center _vals values. */
  // for (int m = 0; m < (int) centers.size(); m++) {
  //   centers[m][0] = centers[m][1] = 0;
  //   for (int n = 0; n < (int) centers_vals[0].size(); n++){
  //     centers_vals[m][n] = 0;
  //   }
  //   center_counts[m] = 0;
  // }
  //
  // // /* Compute the new cluster centers. */
  // for (int l = 0; l < mat.ncol(); l++) {
  //   for (int k = 0; k < mat.nrow(); k++) {
  //     int c_id = clusters[l][k];
  //
  //     if (c_id != -1) {
  //       int ncell = l + (k * mat.ncol());
  //
  //       vector<double> colour;
  //       for (int nval = 0; nval < vals.ncol(); nval++){
  //         double val = vals(ncell, nval);
  //         colour.push_back(val);
  //       }
  //
  //       centers[c_id][0] += k;
  //       centers[c_id][1] += l;
  //       for (int nval = 0; nval < vals.ncol(); nval++){
  //         centers_vals[c_id][nval] += colour[nval];
  //       }
  //
  //       center_counts[c_id] += 1;
  //     }
  //   }
  // }
  //
  // // /* Normalize the clusters. */
  // for (int l = 0; l < (int) centers.size(); l++) {
  //   centers[l][0] /= center_counts[l];
  //   centers[l][1] /= center_counts[l];
  //   for (int nval = 0; nval < vals.ncol(); nval++){
  //     centers_vals[l][nval] /= center_counts[l];
  //   }
  // }
}

writable::integers_matrix Slic::return_clusters(){
  writable::integers_matrix result(clusters[0].size(), clusters.size());

  for (int i = 0; i < clusters.size(); i++) {
    for (int j = 0; j < clusters[0].size(); j++) {
      result(j, i) = clusters[i][j];
    }
  }
  return result;
}

writable::integers_matrix Slic::return_centers(){
  writable::integers_matrix result(centers.size(), 3);
  for (int i = 0; i < (int) centers.size(); i++){
    result(i, 0) = centers[i][0];
    result(i, 1) = centers[i][1];
    result(i, 2) = centers[i][2];
  }
  return result;
}

/*** R
devtools::load_all()
library(landscapemetrics)
library(raster)
library(terra)
library(sf)
library(tmap)

volcanorast = raster(volcano, xmn = 0, xmx = 61, ymn = 0, ymx = 87, crs = 2180)
mode(volcano) = "integer"

b = run_slic(volcano, 5, 5, TRUE, TRUE)
volcanorast2 = raster(b)
extent(volcanorast2) = extent(volcanorast)
vo2 = rast(volcanorast2)
vo3 = terra::as.polygons(vo2, dissolve = TRUE)
vo4 = st_as_sf(vo3, crs = st_crs(volcanorast))
vo5 = st_cast(st_make_valid(vo4), "MULTIPOLYGON")
st_crs(vo5) = st_crs(volcanorast)

tmap_mode("plot")
tm_shape(volcanorast) +
  tm_raster(legend.show = FALSE, style = "cont") +
  tm_shape(vo5, is.master = TRUE) +
  tm_borders()

data("augusta_nlcd")
m = as.matrix(augusta_nlcd)

aaa = run_slic(m, 2, 1, TRUE, TRUE)

augusta_nlcd2 = augusta_nlcd
augusta_nlcd2 = raster(aaa)
extent(augusta_nlcd2) = extent(augusta_nlcd)
an = rast(augusta_nlcd)
an2 = rast(augusta_nlcd2)
an3 = terra::as.polygons(an2, dissolve = TRUE)
an4 = st_as_sf(an3, crs = st_crs(augusta_nlcd))

an5 = st_cast(st_make_valid(an4), "MULTIPOLYGON")
st_crs(an5) = st_crs(augusta_nlcd)

tm_shape(augusta_nlcd) +
  tm_raster(legend.show = FALSE) +
  tm_shape(an5) +
  tm_borders()
*/
