# Script to generate test fixtures with an older version of the package.
library(sf)
library(terra)
library(supercells)
# version 1.0.3
# https://github.com/Nowosad/supercells/tree/285359b8bc862a3fb62d2e3b44ca61772e11d1f6

out_dir = "tests/testthat/testdata"
dir.create("tests/testthat/testdata", showWarnings = FALSE, recursive = TRUE)

v1 = terra::rast(system.file("raster/volcano.tif", package = "supercells"))
v3 = terra::rast(system.file("raster/ortho.tif", package = "supercells"))

v1_results = list(
  v1_a = sf::st_drop_geometry(supercells(v1, 100, compactness = 1)),
  v1_b = sf::st_drop_geometry(supercells(v1, 100, compactness = 1, clean = FALSE)),
  v1_c = sf::st_drop_geometry(supercells(v1, step = 8, compactness = 1)),
  v1_d = sf::st_drop_geometry(supercells(v1, step = 8, compactness = 1, minarea = 3)),
  v1_e = sf::st_drop_geometry(supercells(v1, step = 8, compactness = 0.1, avg_fun = "median", dist_fun = "jaccard")),
  v1_f = sf::st_drop_geometry(supercells(v1, 100, compactness = 1, avg_fun = mean))
)

v3_results = list(
  v3_a = sf::st_drop_geometry(supercells(v3, 100, compactness = 1)),
  v3_b = sf::st_drop_geometry(supercells(v3, 100, compactness = 1, transform = "to_LAB"))
)

saveRDS(v1_results, file.path(out_dir, "v1-supercells-v1.rds"))
saveRDS(v3_results, file.path(out_dir, "v3-supercells-v1.rds"))
