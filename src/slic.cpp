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
  center_counts.clear();
}

void Slic::inits(integers_matrix m){
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
      int colour = m(nrowcenter, ncolcenter);
      vector<int> lm = find_local_minimum(m, nrowcenter, ncolcenter);

      /* Generate the center vector. */
      center.push_back(lm[0]);
      center.push_back(lm[1]);
      center.push_back(colour);

      /* Append to vector of centers. */
      centers.push_back(center);
      center_counts.push_back(0);
    }
  }
}

double Slic::compute_dist(int ci, int y, int x, int value) {

  double dc = sqrt(pow(centers[ci][2] - value, 2));
  double ds = sqrt(pow(centers[ci][0] - x, 2) + pow(centers[ci][1] - y, 2));

  return sqrt(pow(dc / nc, 2) + pow(ds / ns, 2));
}

vector<int> Slic::find_local_minimum(integers_matrix m, int y, int x) {
  int min_grad = -1;

  vector<int> loc_min(2);
  loc_min.at(0) = y;
  loc_min.at(1) = x;

  for (int i = x - 1; i < x + 2; i++) {
    for (int j = y - 1; j < y + 2; j++) {
      int i1 = m(j + 1, i);
      int i2 = m(j, i + 1);
      int i3 = m(j, i);

      /* Compute horizontal and vertical gradients and keep track of the
       minimum. */
      if ((sqrt(pow(i1 - i3, 2)) + sqrt(pow(i2 - i3, 2))) < min_grad) {
        min_grad = fabs(i1 - i3) + fabs(i2 - i3);
        loc_min.at(0) = j;
        loc_min.at(1) = i;
      }
    }
  }
  return loc_min;
}

writable::integers_matrix Slic::generate_superpixels(integers_matrix mat, int step, int nc){
  // cout << "generate_superpixels" << endl;
  this->step = step;
  this->nc = nc;

  /* Clear previous data (if any), and re-initialize it. */
  clear_data();
  inits(mat);

  /* Run EM for 10 iterations (as prescribed by the algorithm). */
  for (int iter = 0; iter < NR_ITERATIONS; iter++) {
    // cout << "iteration: " << NR_ITERATIONS << endl;
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
            int colour = mat(n, m);
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

    /* Clear the center values. */
    for (int m = 0; m < (int) centers.size(); m++) {
      centers[m][0] = centers[m][1] = centers[m][2] = 0;
      center_counts[m] = 0;
    }

    /* Compute the new cluster centers. */
    for (int l = 0; l < mat.ncol(); l++) {
      for (int k = 0; k < mat.nrow(); k++) {
        int c_id = clusters[l][k];

        if (c_id != -1) {
          int colour = mat(k, l);

          centers[c_id][0] += k;
          centers[c_id][1] += l;
          centers[c_id][2] += colour;

          center_counts[c_id] += 1;
        }
      }
    }
    /* Normalize the clusters. */
    for (int l = 0; l < (int) centers.size(); l++) {
      centers[l][0] /= center_counts[l];
      centers[l][1] /= center_counts[l];
      centers[l][2] /= center_counts[l];
    }
  }
  // cout << "centers.size()" << centers.size() << endl;
  create_connectivity(mat);
  clusters = new_clusters;

  /* Clear the center values. */
  for (int m = 0; m < (int) centers.size(); m++) {
    centers[m][0] = centers[m][1] = centers[m][2] = 0;
    center_counts[m] = 0;
  }

  /* Compute the new cluster centers. */
  for (int l = 0; l < mat.ncol(); l++) {
    for (int k = 0; k < mat.nrow(); k++) {
      int c_id = clusters[l][k];

      if (c_id != -1) {
        int colour = mat(k, l);

        centers[c_id][0] += k;
        centers[c_id][1] += l;
        centers[c_id][2] += colour;

        center_counts[c_id] += 1;
      }
    }
  }

  /* Normalize the clusters. */
  for (int l = 0; l < (int) centers.size(); l++) {
    centers[l][0] /= center_counts[l];
    centers[l][1] /= center_counts[l];
    centers[l][2] /= center_counts[l];
  }

  writable::integers_matrix result(centers.size(), 3);
  for (int i = 0; i < (int) centers.size(); i++){
    result(i, 0) = centers[i][0];
    result(i, 1) = centers[i][1];
    result(i, 2) = centers[i][2];
  }
  return result;
}

void Slic::create_connectivity(integers_matrix m) {
  int label = 0;
  int adjlabel = 0;
  const int lims = (m.ncol() * m.nrow()) / ((int)centers.size());
  cout << "lims" << lims << endl;
  const int dx4[4] = {-1,  0,  1,  0};
  const int dy4[4] = { 0, -1,  0,  1};

  for (int i = 0; i < m.ncol(); i++) {
    vector<int> nc;
    for (int j = 0; j < m.nrow(); j++) {
      nc.push_back(-1);
    }
    new_clusters.push_back(nc);
  }

  for (int i = 0; i < m.ncol(); i++) {
    for (int j = 0; j < m.nrow(); j++) {

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

          if (x >= 0 && x < m.ncol() && y >= 0 && y < m.nrow()) {
            if (new_clusters[x][y] >= 0) {
              adjlabel = new_clusters[x][y];
            }
          }
        }

        int count = 1;
        for (int c = 0; c < count; c++) {
          for (int k = 0; k < 4; k++) {
            int x = elements[c][1] + dx4[k], y = elements[c][0] + dy4[k];

            if (x >= 0 && x < m.ncol() && y >= 0 && y < m.nrow()) {
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
      }
    }
  }
}

writable::integers_matrix Slic::colour_with_cluster_means(integers_matrix m) {
  // vector<int> colours(centers.size());

  // /* Gather the colour values per cluster. */
  // for (int i = 0; i < m.ncol(); i++) {
  //   for (int j = 0; j < m.nrow(); j++) {
  //     int index = clusters[i][j];
  //     int colour = m(j, i);
  //
  //     colours[index] = colour;
  //   }
  // }
  //
  // /* Divide by the number of pixels per cluster to get the mean colour. */
  // for (int i = 0; i < (int)colours.size(); i++) {
  //   colours[i] /= center_counts[i];
  // }

  writable::integers_matrix result(m.nrow(), m.ncol());

  /* Fill in. */
  for (int i = 0; i < m.ncol(); i++) {
    for (int j = 0; j < m.nrow(); j++) {
      // int ncolour = colours[clusters[i][j]];
      // result(j, i) = ncolour;
      result(j, i) = clusters[i][j];
      // result(j, i) = new_clusters[i][j];
    }
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
bc = foo(volcano, 5, 5)
bc = st_as_sf(as.data.frame(bc),
              coords = c("V2", "V1"),
              crs = st_crs(volcanorast))


b = foo3(volcano, 5, 5)
volcanorast2 = raster(b)
extent(volcanorast2) = extent(volcanorast)
vo2 = rast(volcanorast2)
vo3 = terra::as.polygons(vo2, dissolve=TRUE)
vo4 = st_as_sf(vo3, crs = st_crs(volcanorast))
vo5 = st_cast(st_make_valid(vo4), "MULTIPOLYGON")
st_crs(vo5) = st_crs(volcanorast)

tmap_mode("plot")
tm_shape(volcanorast) +
  tm_raster(legend.show = FALSE, style = "cont") +
  tm_shape(vo5, is.master = TRUE) +
  tm_borders() +
  tm_shape(bc) +
  tm_dots()

write_sf(vo5, "tmp.gpkg")

data("augusta_nlcd")
m = as.matrix(augusta_nlcd)
a = foo(m, 20, 1)
foo2(m, 1, 1)

aaa = foo3(m, 20, 1)

augusta_nlcd2 = augusta_nlcd
augusta_nlcd2 = raster(aaa)
extent(augusta_nlcd2) = extent(augusta_nlcd)
# plot(augusta_nlcd2)



an = rast(augusta_nlcd)
an2 = rast(augusta_nlcd2)
an3 = terra::as.polygons(an2, dissolve=TRUE)
an4 = st_as_sf(an3, crs = st_crs(augusta_nlcd))

an5 = st_cast(st_make_valid(an4), "MULTIPOLYGON")
st_crs(an5) = st_crs(augusta_nlcd)

library(tmap)
tm_shape(augusta_nlcd) +
  tm_raster(legend.show = FALSE) +
  tm_shape(an5) +
  tm_borders()



library(sf)
library(tmap)
a2 = na.omit(as.data.frame(a))
a2$V1 = a2$V1 * 30 + 1246815
a2$V2 = a2$V2 * 30 + 1249635
a_sf = st_as_sf(a2, coords = c("V2", "V1"), crs = st_crs(augusta_nlcd))

tmap_mode("view")
tm_shape(augusta_nlcd) +
  tm_raster(legend.show = FALSE) +
  tm_shape(a_sf, is.master = TRUE) +
  tm_dots(col = "V3", n = 30)
plot(augusta_nlcd)
starts(200, m)
*/
