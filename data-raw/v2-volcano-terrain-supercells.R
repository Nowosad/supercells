library(terra)
library(sf)
library(supercells)

vol_path = "inst/raster/volcano.tif"
if (!file.exists(vol_path)) {
  vol_path = system.file("raster/volcano.tif", package = "supercells")
}

elevation = terra::rast(vol_path)
slope = terra::terrain(elevation, v = "slope", unit = "radians")
aspect = terra::terrain(elevation, v = "aspect", unit = "radians")
aspect_cos = terra::app(aspect, cos)
aspect_sin = terra::app(aspect, sin)

names(elevation) = "elevation"
names(slope) = "slope"
names(aspect_cos) = "aspect_cos"
names(aspect_sin) = "aspect_sin"

x = c(elevation, slope, aspect_cos, aspect_sin)
vals = terra::values(x, mat = TRUE, na.rm = FALSE)
vals_scaled = scale(vals)
vals_scaled[, "aspect_cos"] = vals_scaled[, "aspect_cos"] / sqrt(2)
vals_scaled[, "aspect_sin"] = vals_scaled[, "aspect_sin"] / sqrt(2)
x_scaled = x
terra::values(x_scaled) = vals_scaled

tune = sc_tune_compactness(x_scaled, step = 8)
regions = sc_slic(x_scaled, step = 8, compactness = tune$compactness)

terra::plot(x_scaled[["elevation"]])
plot(sf::st_geometry(regions), add = TRUE, lwd = 0.2)
