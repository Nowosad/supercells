library(sf)
library(terra)
library(rsi)
terraOptions(tempdir = "~/tmp/", todisk = TRUE)

# aoi = st_transform(
#   st_as_sfc(
#     "POLYGON ((5.670233 52.76554, 5.749025 52.76554, 5.749025 52.731406, 5.670233 52.731406, 5.670233 52.76554))",
#     crs = 4326
#   ),
#   32631
# )

aoi = st_transform(
  st_as_sfc(
    "POLYGON((19.097157 54.1374, 19.097157 54.167358, 19.181614 54.167358, 19.181614 54.1374, 19.097157 54.1374))",
    crs = 4326
  ),
  3035
)

band_map = rsi::sentinel2_band_mapping$planetary_computer_v1[c("B02", "B03", "B04", "B08")]

raw_file = tempfile(fileext = ".tif")
rsi::get_sentinel2_imagery(
  aoi = aoi,
  start_date = "2020-06-24",
  end_date = "2020-06-27",
  rescale_bands = TRUE,
  pixel_x_size = 10,
  pixel_y_size = 10,
  asset_names = band_map,
  composite_function = "median",
  output_filename = raw_file
)

s2 = terra::rast(raw_file)
names(s2) = c("blue", "green", "red", "nir")

terra::writeRaster(
  s2,
  filename = file.path("inst", "raster", "sentinel2_fields.tif"),
  overwrite = TRUE,
  datatype = "FLT4S", # keep float values
  gdal = c(
    "TILED=YES",
    "COMPRESS=ZSTD",
    "ZSTD_LEVEL=22",   # strongest ZSTD
    "PREDICTOR=3",     # best for floating-point rasters
    "NUM_THREADS=ALL_CPUS",
    "BIGTIFF=IF_SAFER"
  )
)
