# remotes::install_github("nowosad/supercells@677715f")
library(supercells)
# devtools::load_all()
library(terra)
library(sf)
vol = rast(system.file("raster/srtm.tif", package = "spDataLarge"))

bench::mark(vol_slic1 = supercells(vol, k = 50, compactness = 1, clean = TRUE),
            vol_slic1 = supercells(vol, k = 50, compactness = 1, clean = TRUE, avg_fun = mean),
            vol_slic1 = supercells(vol, k = 50, compactness = 1, clean = TRUE, avg_fun = "mean2"))

bench::mark(vol_slic1a = supercells(vol, k = 50, compactness = 1, clean = TRUE, avg_fun = "median"),
            vol_slic2a = supercells(vol, k = 50, compactness = 1, clean = TRUE, avg_fun = median),
            check = FALSE)

bench::mark(vol_slic1 = supercells(vol, k = 50, compactness = 1, clean = FALSE),
            vol_slic1 = supercells(vol, k = 50, compactness = 1, clean = FALSE, avg_fun = mean),
            vol_slic1 = supercells(vol, k = 50, compactness = 1, clean = FALSE, avg_fun = "mean2"))

bench::mark(vol_slic1a = supercells(vol, k = 50, compactness = 1, clean = FALSE, avg_fun = "median"),
            vol_slic2a = supercells(vol, k = 50, compactness = 1, clean = FALSE, avg_fun = median),
            check = FALSE)


vol_slic1 = supercells(vol, k = 50, compactness = 1, clean = FALSE)
vol_slic2 = supercells(vol, k = 50, compactness = 1, clean = FALSE, avg_fun = mean)

plot(vol)
plot(st_geometry(vol_slic1), add = TRUE, lwd = 0.2)
plot(st_geometry(vol_slic2), add = TRUE, lwd = 0.2, border = "red")



vol_slic1a = supercells(vol, k = 50, compactness = 1, clean = FALSE, avg_fun = "median")
vol_slic2a = supercells(vol, k = 50, compactness = 1, clean = FALSE, avg_fun = median)

all.equal(vol_slic2a, vol_slic1a)


library(mapview)
mapview(st_geometry(vol_slic1a)) + st_geometry(vol_slic2a)
