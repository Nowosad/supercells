# remotes::install_github("nowosad/supercells@677715f")
library(supercells)
# devtools::load_all()
library(terra)
library(sf)
srtm = rast(system.file("raster/srtm.tif", package = "spDataLarge"))

# benchmarks for different avg_fun specifications
bench::mark(srtm_slic1 = supercells(srtm, k = 50, compactness = 1, clean = TRUE),
            srtm_slic1 = supercells(srtm, k = 50, compactness = 1, clean = TRUE, avg_fun = mean),
            srtm_slic1 = supercells(srtm, k = 50, compactness = 1, clean = TRUE, avg_fun = "mean"))

bench::mark(srtm_slic1a = supercells(srtm, k = 50, compactness = 1, clean = TRUE, avg_fun = "median"),
            srtm_slic2a = supercells(srtm, k = 50, compactness = 1, clean = TRUE, avg_fun = median),
            check = FALSE)

bench::mark(srtm_slic1 = supercells(srtm, k = 50, compactness = 1, clean = FALSE),
            srtm_slic1 = supercells(srtm, k = 50, compactness = 1, clean = FALSE, avg_fun = mean),
            srtm_slic1 = supercells(srtm, k = 50, compactness = 1, clean = FALSE, avg_fun = "mean"))

bench::mark(srtm_slic1a = supercells(srtm, k = 50, compactness = 1, clean = FALSE, avg_fun = "median"),
            srtm_slic2a = supercells(srtm, k = 50, compactness = 1, clean = FALSE, avg_fun = median),
            check = FALSE)

# visual check
srtm_slic1 = supercells(srtm, k = 50, compactness = 1, clean = FALSE)
srtm_slic2 = supercells(srtm, k = 50, compactness = 1, clean = FALSE, avg_fun = mean)

plot(srtm)
plot(st_geometry(srtm_slic1), add = TRUE, lwd = 0.5)
plot(st_geometry(srtm_slic2), add = TRUE, lwd = 0.5, border = "red")

# equality check
srtm_slic1a = supercells(srtm, k = 50, compactness = 1, clean = FALSE, avg_fun = "median")
srtm_slic2a = supercells(srtm, k = 50, compactness = 1, clean = FALSE, avg_fun = median)

all.equal(srtm_slic2a, srtm_slic1a)

