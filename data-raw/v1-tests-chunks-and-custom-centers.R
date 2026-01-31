devtools::load_all()
library(terra)
library(sf)
set.seed(2021-11-05)
ortho = rast(system.file("raster/ortho.tif", package = "supercells"))
y1 = sf::st_sf(geom = sf::st_sample(sf::st_as_sfc(sf::st_bbox(ortho)), 100, type = "random"))

plot(ortho)
plot(st_geometry(y1), add = TRUE)

superpixels0 = supercells(x = ortho,
                         k = y1,
                         step = 100,
                         compactness = 0.5,
                         dist_fun = "jsd",
                         minarea = 16)
plot(superpixels0)

# system.time({superpixels = supercells(x = ortho,
#                                       k = y1,
#                                       step = 50,
#                                       compactness = 0.5,
#                                       dist_fun = "jsd",
#                                       minarea = 16,
#                                       chunks = 150)})
#
