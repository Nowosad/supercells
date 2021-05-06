devtools::load_all()
library(spDataLarge)
library(terra)
library(sf)
library(tmap)
# library(supercell)

get_step = function(r, sp){
  sqrt((ncol(r) * nrow(r)) / sp)
}

# test 1 ------------------------------------------------------------------
srtm = rast(system.file("raster/srtm.tif", package = "spDataLarge"))

dog = raster::raster("data-raw/dog.png"); crs(dog) = "EPSG:2180"
dog2 = rast("data-raw/dog.png"); ext(dog2) = c(0, 320, 0, 240); crs(dog2) = "EPSG:2180"

a = bench::mark(
  srtm_slic = supercell(srtm, 50, 1, "jensen_shannon")
)
a
# A tibble: 1 x 13
# expression      min   median `itr/sec` mem_alloc `gc/sec`
# <bch:expr> <bch:tm> <bch:tm>     <dbl> <bch:byt>    <dbl>
#   1 srtm_slic     733ms    733ms      1.36    26.4MB

b = bench::mark(
  srtm_slic = supercell(dog2, 12, 1, "jensen_shannon"),
  srtm_slic = supercell(dog2, 12, 1, "euclidean"),
  srtm_slic = supercell(dog2, 50, 1, "jensen_shannon"),
  srtm_slic = supercell(dog2, 50, 1, "euclidean"),
  srtm_slic = supercell(dog2, 50, 10, "jensen_shannon"),
  srtm_slic = supercell(dog2, 50, 10, "euclidean"),
  check = FALSE
)
b
# # A tibble: 6 x 13
# expression      min   median `itr/sec` mem_alloc `gc/sec` n_itr
# <bch:expr> <bch:tm> <bch:tm>     <dbl> <bch:byt>    <dbl> <int>
#   1 srtm_slic     9.84s    9.84s     0.102    27.8MB    0.102
#   2 srtm_slic        9s       9s     0.111    27.7MB    0.111
#   3 srtm_slic     6.68s    6.68s     0.150    26.3MB    0.299
#   4 srtm_slic     5.96s    5.96s     0.168    26.4MB    0
#   5 srtm_slic     6.27s    6.27s     0.159    26.2MB    0.159
#   6 srtm_slic     6.07s    6.07s     0.165    26.4MB    0.165
