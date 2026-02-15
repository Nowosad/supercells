# Compare supercells and snic using polygon outputs only.
# Run from the repository root.

if (!requireNamespace("supercells", quietly = TRUE)) {
  stop("Package 'supercells' is required.")
}
if (!requireNamespace("snic", quietly = TRUE)) {
  stop(
    "Package 'snic' is required. Install from GitHub:\n",
    "  remotes::install_github('rolfsimoes/snic')"
  )
}
if (!requireNamespace("terra", quietly = TRUE)) {
  stop("Package 'terra' is required.")
}
if (!requireNamespace("sf", quietly = TRUE)) {
  stop("Package 'sf' is required.")
}
if (!requireNamespace("regional", quietly = TRUE)) {
  stop("Package 'regional' is required for quality comparison metrics.")
}

library(terra)
library(sf)

# Closest cross-package mapping used here:
# - supercells::step  <-> snic::spacing
# - shared compactness numeric value
# - supercells clean = FALSE to avoid extra post-segmentation cleaning
# - snic seeds: rectangular grid with spacing/padding from README examples
# - regional::reg_inhomogeneity/reg_isolation used for quality comparison

make_snic_seeds <- function(img, spacing, padding = floor(spacing / 2)) {
  snic::snic_grid(
    x = img,
    type = "rectangular",
    spacing = spacing,
    padding = padding
  )
}

run_pair <- function(img, spacing, compactness, clean = FALSE, iter = 10) {
  sc_poly <- supercells::sc_slic(
    x = img,
    step = spacing,
    compactness = compactness,
    clean = clean,
    iter = iter
  )
  seeds <- make_snic_seeds(img, spacing = spacing)
  snic_obj <- snic::snic(
    x = img,
    seeds = seeds,
    compactness = compactness
  )

  snic_ids <- snic::snic_get_seg(snic_obj)
  snic_poly <- sf::st_as_sf(
    terra::as.polygons(snic_ids, dissolve = TRUE, na.rm = TRUE)
  )

  list(
    sc_poly = sc_poly,
    snic_poly = snic_poly
  )
}

add_regional_metrics <- function(polygons, raster, sample_size = 1) {
  polygons$inh <- regional::reg_inhomogeneity(
    region = polygons,
    raster = raster,
    sample_size = sample_size
  )
  polygons$iso <- regional::reg_isolation(
    region = polygons,
    raster = raster,
    sample_size = sample_size
  )
  polygons
}

quality_summary <- function(name, polygons) {
  data.frame(
    method = name,
    n = nrow(polygons),
    inhomogeneity_mean = mean(polygons$inh, na.rm = TRUE),
    inhomogeneity_median = stats::median(polygons$inh, na.rm = TRUE),
    isolation_mean = mean(polygons$iso, na.rm = TRUE),
    isolation_median = stats::median(polygons$iso, na.rm = TRUE)
  )
}

cat("\n--- Example 1: single-band volcano raster ---\n")
vol <- terra::rast("inst/raster/volcano.tif")
ex1 <- run_pair(vol, spacing = 8, compactness = 0.5, clean = TRUE, iter = 10)
cat("supercells segments:", nrow(ex1$sc_poly), "\n")
cat("snic segments:", nrow(ex1$snic_poly), "\n")

op <- par(mfrow = c(1, 2), mar = c(2, 2, 2, 2))
terra::plot(vol, main = "supercells polygons (volcano)")
plot(sf::st_geometry(ex1$sc_poly), add = TRUE, border = "red", lwd = 0.7)
terra::plot(vol, main = "snic polygons (volcano)")
plot(sf::st_geometry(ex1$snic_poly), add = TRUE, border = "yellow", lwd = 0.7)
par(op)

ex1$sc_poly <- add_regional_metrics(ex1$sc_poly, vol, sample_size = 1)
ex1$snic_poly <- add_regional_metrics(ex1$snic_poly, vol, sample_size = 1)
ex1_quality <- rbind(
  quality_summary("supercells", ex1$sc_poly),
  quality_summary("snic", ex1$snic_poly)
)
cat("\nQuality summary (volcano):\n")
print(ex1_quality, row.names = FALSE)

cat("\n--- Example 2: Sentinel-2 toy (multi-band, RGB) ---\n")
s2 <- terra::rast("data-raw/sentinel2_fields_raw_multipolygon.tif")[[1:3]]
ex2 <- run_pair(s2, spacing = 20, compactness = 0.5, clean = TRUE, iter = 10)
cat("supercells segments:", nrow(ex2$sc_poly), "\n")
cat("snic segments:", nrow(ex2$snic_poly), "\n")

op <- par(mfrow = c(1, 2), mar = c(2, 2, 2, 2))
terra::plotRGB(s2, r = 3, g = 2, b = 1, stretch = "lin", main = "supercells (RGB)")
plot(sf::st_geometry(ex2$sc_poly), add = TRUE, border = "red", lwd = 0.7)
terra::plotRGB(s2, r = 3, g = 2, b = 1, stretch = "lin", main = "snic polygons (RGB)")
plot(sf::st_geometry(ex2$snic_poly), add = TRUE, border = "yellow", lwd = 0.7)
par(op)

ex2$sc_poly <- add_regional_metrics(ex2$sc_poly, s2, sample_size = 0.2)
ex2$snic_poly <- add_regional_metrics(ex2$snic_poly, s2, sample_size = 0.2)
ex2_quality <- rbind(
  quality_summary("supercells", ex2$sc_poly),
  quality_summary("snic", ex2$snic_poly)
)
cat("\nQuality summary (sentinel2 RGB, sample_size = 0.2):\n")
print(ex2_quality, row.names = FALSE)
