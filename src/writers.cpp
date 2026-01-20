#include "slic_core.h"
#include "cpp11.hpp"
#include "cpp11/integers.hpp"
#include "cpp11/matrix.hpp"

cpp11::writable::integers_matrix<> return_clusters(const SlicCore& slic) {
  const auto& clusters = slic.clusters_ref();
  int isize = clusters.size();
  int jsize = clusters[0].size();

  cpp11::writable::integers_matrix<> result(jsize, isize);

  for (int i = 0; i < isize; i++) {
    for (int j = 0; j < jsize; j++) {
      result(j, i) = clusters[i][j];
    }
  }
  return result;
}

cpp11::writable::doubles_matrix<> return_centers(const SlicCore& slic) {
  const auto& centers = slic.centers_ref();
  cpp11::writable::doubles_matrix<> result(centers.size(), 2);
  for (int i = 0; i < (int) centers.size(); i++) {
    result(i, 1) = centers[i][0]; /*y*/
    result(i, 0) = centers[i][1]; /*x*/
  }
  return result;
}

cpp11::writable::doubles_matrix<> return_centers_vals(const SlicCore& slic) {
  const auto& centers_vals = slic.centers_vals_ref();
  const auto& mat_dims = slic.mat_dims_ref();
  cpp11::writable::doubles_matrix<> result(centers_vals.size(), mat_dims[2]);
  for (int i = 0; i < (int) centers_vals.size(); i++) {
    for (int nval = 0; nval < mat_dims[2]; nval++) {
      result(i, nval) = centers_vals[i][nval];
    }
  }
  return result;
}
