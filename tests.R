devtools::load_all()
library(spDataLarge)
library(terra)
library(sf)
library(tmap)
# library(slicr)

# test 1 ------------------------------------------------------------------
srtm = rast(system.file("raster/srtm.tif", package = "spDataLarge"))
srtmm = sqrt(as.matrix(srtm, wide = TRUE))
mode(srtmm) = "integer"
srtm_slic = run_slic(srtmm, 50, 5, FALSE, TRUE)
srtm_slic = rast(srtm_slic)
ext(srtm_slic) = ext(srtm)
crs(srtm_slic) = st_crs(srtm)$proj4string
srtm_slic = st_as_sf(terra::as.polygons(srtm_slic, dissolve = TRUE))

tm_shape(srtm) +
  tm_raster(style = "cont") +
  tm_shape(srtm_slic) +
  tm_borders()

# test 2 ------------------------------------------------------------------
nlcd = rast(system.file("raster/nlcd2011.tif", package = "spDataLarge"))
nlcdm = as.matrix(nlcd, wide = TRUE)
mode(nlcdm) = "integer"
nlcd_slic = run_slic(nlcdm, 500, 5, TRUE, TRUE)
nlcd_slic = rast(nlcd_slic)
ext(nlcd_slic) = ext(nlcd)
crs(nlcd_slic) = st_crs(nlcd)$proj4string
nlcd_slic = st_as_sf(terra::as.polygons(nlcd_slic, dissolve = TRUE))

tm_shape(nlcd) +
  tm_raster(legend.show = FALSE, style = "cat") +
  tm_shape(nlcd_slic) +
  tm_borders()

# test 3 ------------------------------------------------------------------
landsat = rast(system.file("raster/landsat.tif", package = "spDataLarge"))
plotRGB(landsat, r = 3, g = 2, b = 1, stretch = "lin")

run_slic()
