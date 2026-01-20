library(terra)
library(sf)
library(supercells)
library(bench)

library(rsi)
my_point = st_sfc(st_point(c(22.850000, 54.250000)), crs = "EPSG:4326")
my_point = st_transform(my_point, crs = "EPSG:2180")
my_poly = st_sf(geom = st_buffer(my_point, dist = 5000))
sa_landsat = get_landsat_imagery(my_poly,
                                 start_date = "2023-09-01", end_date = "2023-10-31",
                                 output_filename = tempfile(fileext = ".tif"))
lnd_low = rast(sa_landsat)
lnd_mid = disagg(lnd_low, 4)
lnd_high = disagg(lnd_mid, 8)

results_high = bench::mark(
        s = supercells(lnd_mid, compactness = 1, step = 30, clean = FALSE, iter = 1, diagnostics = TRUE),
        iterations = 10
)

s = supercells(lnd_mid, compactness = 1, step = 30, clean = FALSE, iter = 10, diagnostics = TRUE)
plot(s)
a = attr(s, "diagnostics")
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

df_cluster = as_tibble(a$per_cluster)
sc = cbind(s, na.omit(df_cluster))
tm_shape(sc) +
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

v = rast(lnd_mid, nlyrs = length(pp))

for (i in seq_along(pp)) {
  values(v[[i]]) = as.vector(t(pp[[i]]))
}

plot(v[[1]])
plot(s[0], add = TRUE, lwd = 1, border = "red")
plot(v[[2]])
plot(s[0], add = TRUE, lwd = 1, border = "red")
plot(v[[3]])
plot(s[0], add = TRUE, lwd = 1, border = "red")

