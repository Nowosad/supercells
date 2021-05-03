devtools::load_all()
library(spDataLarge)
library(terra)
library(sf)
library(tmap)
# library(slicr)

get_step = function(r, sp){
  sqrt((ncol(r) * nrow(r)) / sp)
}

# test 1 ------------------------------------------------------------------
srtm = rast(system.file("raster/srtm.tif", package = "spDataLarge"))
srtmm = as.matrix(srtm, wide = TRUE)
srtmv = as.matrix(as.data.frame(srtm, cell = TRUE)[2])
mode(srtmm) = "integer"
mode(srtmv) = "integer"
srtm_slic = run_slic(srtmm, srtmv, 23.04913, 50, FALSE, TRUE)
# srtm_slic = slicr2:::run_slic(srtmm, 50, 5, TRUE, TRUE)

srtm_slic = rast(srtm_slic)
ext(srtm_slic) = ext(srtm)
crs(srtm_slic) = st_crs(srtm)$proj4string
srtm_slic = st_as_sf(terra::as.polygons(srtm_slic, dissolve = TRUE))

# a = st_collection_extract(srtm_slic, "POLYGON")
# b = st_cast(a, "POLYGON")

tm_shape(srtm) +
  tm_raster(style = "cont") +
  tm_shape(srtm_slic) +
  tm_borders()

# test 2 ------------------------------------------------------------------
nlcd = rast(system.file("raster/nlcd2011.tif", package = "spDataLarge"))
nlcdm = as.matrix(nlcd, wide = TRUE)
mode(nlcdm) = "integer"
nlcd_slic = run_slic(nlcdm, 150, 100, TRUE, TRUE)
nlcd_slic = rast(nlcd_slic)
ext(nlcd_slic) = ext(nlcd)
crs(nlcd_slic) = st_crs(nlcd)$proj4string
nlcd_slic = st_as_sf(terra::as.polygons(nlcd_slic, dissolve = TRUE))

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

plotRGB(landsat2, r = 3, g = 2, b = 1)

landsatm = as.matrix(landsat2[[1]], wide = TRUE)
landsatv = as.matrix(as.data.frame(landsat2, cell = TRUE)[-c(1, 5)])
mode(landsatm) = "integer"
mode(landsatv) = "integer"
landsat_slic = run_slic(landsatm, landsatv, 50, 50, TRUE, TRUE)
landsat_slic = rast(landsat_slic)
ext(landsat_slic) = ext(landsat)
crs(landsat_slic) = st_crs(landsat)$proj4string
landsat_slict = terra::as.polygons(landsat_slic, dissolve = TRUE)
landsat_slic = st_as_sf(landsat_slict)


plotRGB(landsat2, r = 3, g = 2, b = 1)
lines(landsat_slict, add = TRUE, lwd = 3, col = "red")

# terra::writeRaster(landsat2, "l.tif")
# write_sf(landsat_slic, "l.gpkg")
# plot(landsat_slic)

l1 = dplyr::filter(stars::st_as_stars(landsat), band == 3)
tm_shape(l1) +
  tm_raster(style = "cont") +
  tm_shape(landsat_slic) +
  tm_borders()


# test 4 ------------------------------------------------------------------
dog = rast("dog.png")
ext(dog) = c(0, 240, 0, 320)
plotRGB(dog, r = 3, g = 2, b = 1, stretch = "lin")

dogm = as.matrix(dog[[1]], wide = TRUE)
dogv = as.matrix(as.data.frame(dog, cell = TRUE)[-c(1, 5)])
mode(dogm) = "integer"
mode(dogv) = "integer"
dog_slic = run_slic(dogm, dogv, 13.85, 200, TRUE, TRUE)
dog_slic = rast(dog_slic)
ext(dog_slic) = ext(dog)
crs(dog_slic) = st_crs(dog)$proj4string
dog_slict = terra::as.polygons(dog_slic, dissolve = TRUE)
dog_slic = st_as_sf(dog_slict)

plotRGB(dog, r = 3, g = 2, b = 1, stretch = "lin")
lines(dog_slict, add = TRUE, lwd = 3, col = "red")
# plot(dog_slic)
