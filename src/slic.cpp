#include "slic.hpp"

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
  // cout << "distances.size()" << distances.size() << endl;
  // cout << "clusters.size()" << clusters.size() << endl;
  // cout << "step" << step << endl;
  // cout << "m.ncol() - step/2)" << m.ncol() - step/2 << endl;

  for (int ncolcenter = step; ncolcenter < m.ncol() - step/2; ncolcenter += step){
    for (int nrowcenter = step; nrowcenter < m.nrow() - step/2; nrowcenter += step){
      vector<double> center;
      int colour = m(nrowcenter, ncolcenter);
      cout << "i" << ncolcenter << endl;
      cout << "j" << nrowcenter << endl;

      vector<int> lm = find_local_minimum(m, nrowcenter, ncolcenter);

      // cout << "colour" << colour << endl;
      // cout << "lm0" << lm[0][0] << endl;
      // cout << "lm1" << lm[0][1] << endl;

      /* Generate the center vector. */
      center.push_back(lm[0]);
      center.push_back(lm[1]);
      center.push_back(colour);

      /* Append to vector of centers. */
      centers.push_back(center);
      center_counts.push_back(0);
    }
  }
  // cout << "centers.size()" << centers.size() << endl;
  // cout << "center_counts.size()" << center_counts.size() << endl;
  //
  // cout << "inits ends" << endl;
}

double Slic::compute_dist(int ci, int y, int x, int value) {

  double dc = sqrt(pow(centers[ci][2] - value, 2));
  double ds = sqrt(pow(centers[ci][0] - x, 2) + pow(centers[ci][1] - y, 2));

  return sqrt(pow(dc / nc, 2) + pow(ds / ns, 2));
}

vector<int> Slic::find_local_minimum(integers_matrix m, int y, int x) {
  int min_grad = -1;
  // cout << "min_grad" << min_grad << endl;

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
  cout << "generate_superpixels" << endl;
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
      // centers[l][3] /= center_counts[l];
      // centers[l][4] /= center_counts[l];
    }
  }
  // cout << "centers.size()" << centers.size() << endl;

  writable::integers_matrix result(centers.size(), 3);
  for (int i = 0; i < (int) centers.size(); i++){
    result(i, 0) = centers[i][0];
    result(i, 1) = centers[i][1];
    result(i, 2) = centers[i][2];
  }
  return result;
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
    }
  }
  return result;
}

// void starts(int step, integers_matrix x){
//   int x_ncol = x.ncol();
//   int x_nrow = x.nrow();
//   int x_ncol_step = step / 2;
//   int output_rows =  floor(x_ncol / step) * floor(x_nrow / step);
//   integers_matrix =
//   for (x_ncol_step; x_ncol_step < x_ncol; x_ncol_step += step){
//     cout << x_ncol_step << "- xcol" << endl;
//     int x_nrow_step = step / 2;
//     for (x_nrow_step; x_nrow_step < x_nrow; x_nrow_step += step){
//       cout << x_nrow_step << "- xrow"<< endl;
//     }
//   }
//   // integers_matrix r(0,0);
// }


  // def initial_cluster_center(S,img,img_h,img_w,clusters):
  //   h = S // 2
  //   w = S // 2
  //   while h < img_h:
  //     while w < img_w:
  //       clusters.append(make_superPixel(h, w,img))
  //       w += S
  //       w = S // 2
  //       h += S
  //       return clusters



/*** R
library(landscapemetrics)
library(raster)
data("augusta_nlcd")
m = as.matrix(augusta_nlcd)
a = foo(m, 20, 1)
foo2(m, 1, 1)

aaa = foo3(m, 20, 1)

augusta_nlcd2 = augusta_nlcd
augusta_nlcd2 = raster(aaa)
extent(augusta_nlcd2) = extent(augusta_nlcd)
plot(augusta_nlcd2)

library(terra)
library(sf)
library
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

volcanorast = raster(volcano, xmn = 0, xmx = 61, ymn = 0, ymx = 87)
mode(volcano) = "integer"
b = foo3(volcano, 5, 1)
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
