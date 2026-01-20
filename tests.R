devtools::load_all()
library(terra)
library(tidyr)
library(ggplot2)
library(tmap)
# One variable

vol = terra::rast(system.file("raster/volcano.tif", package = "supercells"))
vol2 = disagg(vol, 2)
system.time({vol_slic1 = supercells(vol2, k = 10, compactness = 46.5, diagnostics = FALSE, clean = TRUE)})
system.time({vol_slic1 = supercells(vol2, k = 100, compactness = 1, diagnostics = TRUE, clean = TRUE)})
terra::plot(vol)
plot(sf::st_geometry(vol_slic1), add = TRUE, lwd = 0.2)

a = attr(vol_slic1, "diagnostics")
str(a)

# Iteration metrics
df = as_tibble(a$iteration)
df$iter = 1:nrow(df)
df_long = pivot_longer(df, cols = -iter, names_to = "metric", values_to = "value")
ggplot(df_long, aes(x = as.integer(iter), y = value, color = metric)) +
  geom_line() +
  facet_wrap(metric~., scales = "free_y", ncol = 1) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
  theme_minimal()

# Global metrics
# The remaining scale (global) statistics are the means of the two distance components, before and after weighting:

# spatial_mean: average raw spatial distance from each pixel to its assigned center (in pixel units).
# value_mean: average raw value‑space distance from each pixel to its assigned center (in the units of your distance function).
# weighted_spatial_mean: average spatial distance after applying the spatial weight (spatial/step).
# weighted_value_mean: average value distance after applying the value weight (value/compactness).
# Interpretation: if weighted_value_mean >> weighted_spatial_mean, the value term dominates assignments; if it’s smaller, spatial proximity dominates.
gdf = a$scale
gdf

# Cluster metrics
# per_cluster: summaries across pixels in each supercell.

# mean_value: value‑space distance to center.
# mean_spatial: spatial distance to center.
# ratio_mean: (value/compactness) / (spatial/step) using mean; >1 means value term dominates.

vol_slic1c = cbind(vol_slic1, a$per_cluster)
tm_shape(vol_slic1c) +
  tm_polygons(c("mean_value", "mean_spatial", "ratio_mean"),
              fill.scale = tm_scale_continuous(values = "viridis"),
              fill.legend = tm_legend(title = "")) +
  tm_facets(nrow = 1)

# Per pixel
# Per‑pixel diagnostics are three matrices with the same row/column layout as the original raster (rows = y, cols = x). Each entry corresponds to that cell’s assigned center:

# per_pixel$spatial: Euclidean pixel distance from the cell to its cluster center (in pixel units).
# per_pixel$value: Distance in feature space (using your selected dist_fun) from the cell’s values to the center’s average values.
# per_pixel$combined: The weighted distance used for assignment: sqrt((value/compactness)^2 + (spatial/step)^2); NA where the cell is NA or unassigned.
# These are useful for mapping where assignments are driven by space vs values (e.g., high value but low spatial indicates a value outlier close to the center).
pp = a$per_pixel

v = rast(vol2, nlyrs = length(pp))

for (i in seq_along(pp)) {
  values(v[[i]]) = as.vector(t(pp[[i]]))
}

plot(v[[1]])
plot(vol_slic1[0], add = TRUE, lwd = 1, border = "red")
plot(v[[2]])
plot(vol_slic1[0], add = TRUE, lwd = 1, border = "red")
plot(v[[3]])
plot(vol_slic1[0], add = TRUE, lwd = 1, border = "red")
