devtools::load_all()
library(terra)
library(sf)
# library(supercells)

ortho = rast(system.file("raster/ortho.tif", package = "supercells"))
ortho2 = extend(ortho, c(2000, 2000))

system.time({superpixels = supercells(x = ortho2,
                                      step = 20,
                                      compactness = 1,
                                      minarea = 16,
                                      chunks = 250)})

library(future)
plan(multisession, workers = 3)

system.time({superpixels2 = supercells(x = ortho2,
                                      step = 20,
                                      compactness = 1,
                                      minarea = 16,
                                      chunks = 250,
                                      future = TRUE)})

system.time({superpixels2 = supercells(x = ortho2,
                                       step = 20,
                                       compactness = 0.5,
                                       minarea = 16,
                                       chunks = 500,
                                       future = TRUE)})

library(future)
plan(multisession, workers = 3)
#
system.time({superpixels2 = supercells(x = ortho2,
                                       step = 20,
                                       compactness = 0.5,
                                       minarea = 16,
                                       chunks = 500,
                                       future = FALSE)})

system.time({superpixels2 = supercells(x = ortho2,
                                       step = 20,
                                       compactness = 0.5,
                                       minarea = 16,
                                       chunks = 500,
                                       future = TRUE)})
