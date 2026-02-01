devtools::load_all()
library(terra)
library(sf)

v1 = rast(system.file("raster/volcano.tif", package = "supercells"))
v3 = rast(system.file("raster/ortho.tif", package = "supercells"))

plot_merge = function(r, sc, merged, title, rgb = FALSE) {
  if (isTRUE(rgb)) {
    terra::plotRGB(r, main = title)
  } else {
    plot(r, main = title)
  }
  # plot(st_geometry(sc), add = TRUE, lwd = 0.5)
  plot(st_geometry(merged), add = TRUE, lwd = 2, border = "red")
}

cat("v1: base segmentation\n")
v1_sc = sc_slic(v1, step = 8, compactness = 1, metadata = TRUE, clean = TRUE)
cat("  n =", nrow(v1_sc), "\n")

cat("v1: greedy merges (target_k)\n")
v1_greedy_k = sc_merge_supercells(v1_sc, dist_fun = "euclidean", method = "greedy",
                                  method_opts = list(target_k = 25),
                                  weight = "area")
cat("  n =", nrow(v1_greedy_k), "\n")
plot_merge(v1, v1_sc, v1_greedy_k, "v1: greedy (target_k)")

cat("v1: greedy merges (tau)\n")
v1_greedy_tau = sc_merge_supercells(v1_sc, dist_fun = "euclidean", method = "greedy",
                                    method_opts = list(tau = 5.5),
                                    weight = "area")
cat("  n =", nrow(v1_greedy_tau), "\n")
plot_merge(v1, v1_sc, v1_greedy_tau, "v1: greedy (tau)")

cat("v1: FH merges (kappa)\n")
v1_fh = sc_merge_supercells(v1_sc, dist_fun = "euclidean", method = "fh",
                            method_opts = list(kappa = 6))
cat("  n =", nrow(v1_fh), "\n")
plot_merge(v1, v1_sc, v1_fh, "v1: FH (kappa)")

cat("v1: MST merges (target_k)\n")
v1_mst_k = sc_merge_supercells(v1_sc, dist_fun = "euclidean", method = "mst",
                               method_opts = list(target_k = 25),
                               weight = "area")
cat("  n =", nrow(v1_mst_k), "\n")
plot_merge(v1, v1_sc, v1_mst_k, "v1: MST (target_k)")

cat("v3: base segmentation (RGB to LAB)\n")
v3_vals = terra::values(v3, mat = TRUE, na.rm = FALSE)
v3_vals = v3_vals / 255
v3_lab = grDevices::convertColor(v3_vals, from = "sRGB", to = "Lab")
v3_lab_r = v3
terra::values(v3_lab_r) = v3_lab

v3_sc = sc_slic(v3_lab_r, step = 15, compactness = 20, metadata = TRUE, clean = TRUE)
cat("  n =", nrow(v3_sc), "\n")

cat("v3: greedy merges (target_k)\n")
v3_greedy_k = sc_merge_supercells(v3_sc, dist_fun = "euclidean", method = "greedy",
                                  method_opts = list(target_k = 78),
                                  weight = "area")
cat("  n =", nrow(v3_greedy_k), "\n")
plot_merge(v3[[1:3]], v3_sc, v3_greedy_k, "v3: greedy (target_k)", rgb = TRUE)

cat("v3: greedy merges (tau)\n")
v3_greedy_tau = sc_merge_supercells(v3_sc, dist_fun = "euclidean", method = "greedy",
                                    method_opts = list(tau = 8),
                                    weight = "area")
cat("  n =", nrow(v3_greedy_tau), "\n")
plot_merge(v3[[1:3]], v3_sc, v3_greedy_tau, "v3: greedy (tau)", rgb = TRUE)

cat("v3: FH merges (kappa)\n")
v3_fh = sc_merge_supercells(v3_sc, dist_fun = "euclidean", method = "fh",
                            method_opts = list(kappa = 17.7))
cat("  n =", nrow(v3_fh), "\n")
plot_merge(v3[[1:3]], v3_sc, v3_fh, "v3: FH (kappa)", rgb = TRUE)

cat("v3: MST merges (target_k)\n")
v3_mst_k = sc_merge_supercells(v3_sc, dist_fun = "euclidean", method = "mst",
                               method_opts = list(target_k = 78),
                               weight = "area")
cat("  n =", nrow(v3_mst_k), "\n")
plot_merge(v3[[1:3]], v3_sc, v3_mst_k, "v3: MST (target_k)", rgb = TRUE)


bench::mark(
  v3_greedy_k = sc_merge_supercells(v3_sc, dist_fun = "euclidean", method = "greedy",
                                    method_opts = list(target_k = 78), weight = "area"),
  v3_greedy_tau = sc_merge_supercells(v3_sc, dist_fun = "euclidean", method = "greedy",
                                     method_opts = list(tau = 8), weight = "area"),
  v3_fh = sc_merge_supercells(v3_sc, dist_fun = "euclidean", method = "fh",
                              method_opts = list(kappa = 17.7)),
  v3_mst_k = sc_merge_supercells(v3_sc, dist_fun = "euclidean", method = "mst",
                                 method_opts = list(target_k = 78), weight = "area"),
  check = FALSE
)
