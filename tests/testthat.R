library(testthat)
library(supercells)
library(terra)

v1 = rast(system.file("raster/volcano.tif", package = "supercells"))
v3 = rast(system.file("raster/ortho.tif", package = "supercells"))


