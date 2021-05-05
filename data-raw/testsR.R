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
srtm_slic = supercell(srtm, 23.04913, 1, "jensen_shannon")

tm_shape(srtm) +
  tm_raster(style = "cont") +
  tm_shape(srtm_slic) +
  tm_borders()

# test 2 ------------------------------------------------------------------
nlcd = rast(system.file("raster/nlcd2011.tif", package = "spDataLarge"))
nlcd_slic = supercell(nlcd, 150, 100)

tm_shape(nlcd) +
  tm_raster(legend.show = FALSE, drop.levels = TRUE) +
  tm_shape(nlcd_slic) +
  tm_borders(col = "black")

# test 3 ------------------------------------------------------------------
landsat = rast(system.file("raster/landsat.tif", package = "spDataLarge"))

.linStretch <- function (x) {
  v_vals = terra::values(x)
  v <- stats::quantile(v_vals, c(0.02, 0.98), na.rm = TRUE)
  temp <- (255 * (x - v[1]))/(v[2] - v[1])
  temp[temp < 0] <- 0
  temp[temp > 255] <- 255
  return(temp)
}
landsat2 = landsat
landsat2[[1]] = .linStretch(landsat[[1]])
landsat2[[2]] = .linStretch(landsat[[2]])
landsat2[[3]] = .linStretch(landsat[[3]])
landsat2[[4]] = .linStretch(landsat[[4]])

landsat_slic = supercell(landsat2, 150, 10, "jensen_shannon")

tm_shape(landsat2) +
  tm_rgb(r = 3, g = 2, b = 1, legend.show = FALSE, drop.levels = TRUE) +
  tm_shape(landsat_slic) +
  tm_borders(col = "red")

# test 4 ------------------------------------------------------------------
dog = raster::raster("dog.png"); crs(dog) = "EPSG:2180"
dog2 = rast("dog.png"); ext(dog2) = c(0, 320, 0, 240); crs(dog2) = "EPSG:2180"
dog_slic = supercell(dog2, 12, 100)

terra::plotRGB(dog2, r = 3, g = 2, b = 1, stretch = "lin")
lines(vect(dog_slic), add = TRUE, col = "red")

tm_shape(dog) +
  tm_rgb(r = 3, g = 2, b = 1, legend.show = FALSE, drop.levels = TRUE) +
  tm_shape(dog_slic) +
  tm_borders(col = "red")
