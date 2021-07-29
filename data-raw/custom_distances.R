devtools::load_all()
library(spDataLarge)
library(terra)
library(sf)
library(tmap)

# test 1 ------------------------------------------------------------------
srtm = rast(system.file("raster/srtm.tif", package = "spDataLarge"))

srtm_slic1 = supercells(srtm, 50, 0.01, "euclidean", iter = 10)
system.time({srtm_slic2 = supercells(srtm, 50, 0.0000000001, "jensonshannon", iter = 10)})
system.time({srtm_slic3 = supercells(srtm, 50, 0.0000000001, "jenson_shannon", iter = 10)})
srtm_slic4 = supercells(srtm, 50, 0.01, "jaccard", iter = 10)
