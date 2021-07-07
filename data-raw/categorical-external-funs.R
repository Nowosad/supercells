devtools::load_all()
library(motif)
library(sf)
library(terra)


# small test
# to make sense of it you need to uncomment printing in slic.cpp
r = rast(matrix(c(rep(1, 15), 10, rep(2, 9)), nrow = 5), crs = "EPSG:2180")
plot(r)

vol_slic1a = supercells(r, step = 2, compactness = 1, clean = FALSE, avg_fun = "median", iter = 30)

plot(t(r))
plot(st_geometry(vol_slic1a), add = TRUE, lwd = 0.5)

# 1 -----------------------------------------------------------------------
nlcd = rast(system.file("raster/nlcd2011.tif", package = "spDataLarge"))

bench::mark(vol_slic1a = supercells(nlcd, k = 50, compactness = 1, clean = TRUE, avg_fun = "median"),
            vol_slic2a = supercells(nlcd, k = 50, compactness = 1, clean = TRUE, avg_fun = median),
            check = FALSE)


vol_slic1a = supercells(nlcd, k = 50, compactness = 1, clean = TRUE, avg_fun = "median")
vol_slic2a = supercells(nlcd, k = 50, compactness = 1, clean = TRUE, avg_fun = "mean")

plot(nlcd)
plot(st_geometry(vol_slic1a), add = TRUE, lwd = 0.5)
plot(st_geometry(vol_slic2a), add = TRUE, lwd = 0.5, border = "red")
