#include "slic.h"

writable::integers_matrix<> Slic::return_clusters(){
  int isize = clusters.size();
  int jsize = clusters[0].size();

  writable::integers_matrix<> result(jsize, isize);

  for (int i = 0; i < isize; i++) {
    for (int j = 0; j < jsize; j++) {
      result(j, i) = clusters[i][j];
    }
  }
  return result;
}

writable::doubles_matrix<> Slic::return_centers(){
  writable::doubles_matrix<> result(centers.size(), 2);
  for (int i = 0; i < (int) centers.size(); i++){
    result(i, 1) = centers[i][0]; /*y*/
    result(i, 0) = centers[i][1]; /*x*/
  }
  return result;
}

writable::doubles_matrix<> Slic::return_centers_vals(){
  writable::doubles_matrix<> result(centers_vals.size(), mat_dims[2]);
  for (int i = 0; i < (int) centers_vals.size(); i++){
    for (int nval = 0; nval < mat_dims[2]; nval++){
      result(i, nval) = centers_vals[i][nval];
    }
  }
  return result;
}