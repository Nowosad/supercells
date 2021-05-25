devtools::load_all()
library(spDataLarge)
library(terra)
library(sf)
library(stars)
library(ggspatial)
library(ggplot2)
library(gganimate)
library(colorspace)
library(purrr)
# library(supercells)

# test 1 ------------------------------------------------------------------
srtm = rast(system.file("raster/srtm.tif", package = "spDataLarge"))
srtm_slic = supercells(srtm, 100, 1, "jensen_shannon")
ggplot() +
  geom_stars(data = stars::st_as_stars(raster::raster(srtm))) +
  scale_fill_continuous_sequential("Terrain 2") +
  geom_sf(data = srtm_slic, fill = NA, color = "red", size = 0.1) +
  stars:::theme_stars() +
  theme(legend.position = "none") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = NULL, y = NULL)


# euclidean ---------------------------------------------------------------
supercell2 = function(step, x, nc, dist_fun = "euclidean", clean = TRUE){
  supercells(x = x, k = step, compactness = nc, dist_fun = dist_fun, clean = clean)
}
add_cols = function(x, col_value){
  x$sizes = as.factor(col_value)
  x
}
supercell3 = function(nc, x, step, dist_fun = "euclidean", clean = TRUE, transform = NULL){
  supercells(x = x, k = step, compactness = nc, dist_fun = dist_fun, clean = clean, transform = transform)
}

# euclidean ---------------------------------------------------------------
nc = c(seq(400, 100, by = -15), seq(90, 10, by = -5), seq(9, 1, by = -1))
names(nc) = nc
system.time({srtm_slic3 = map(nc, supercell3, srtm, 200, "euclidean")})
srtm_slic3b = map2(srtm_slic3, nc, add_cols)
srtm_slic3c = do.call(rbind, srtm_slic3b)

g3a = ggplot() +
  geom_stars(data = stars::st_as_stars(raster::raster(srtm))) +
  scale_fill_continuous_sequential("Terrain 2") +
  geom_sf(data = srtm_slic3c, fill = NA, color = "#DA3DBE", size = 0.1, aes(group = sizes)) +
  stars:::theme_stars() +
  transition_states(sizes, state_length = 2) +
  ease_aes("bounce-in") +
  enter_fade() +
  exit_shrink() +
  theme(legend.position = "none") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = NULL, y = NULL, title = "Compactness: {closest_state}")

# animate(g1a, renderer = ffmpeg_renderer())
anim_save("g3a2.gif", g3a, height = 537/100, width = 420/100, units = "in", res = 150)


# rgb ---------------------------------------------------------------------
v3 = rast(system.file("raster/ortho.tif", package = "supercells"))

nc = c(seq(100, 50, by = -10), seq(50, 10, by = -2), seq(9, 1, by = -1), 0.5, 0.1)
names(nc) = nc
system.time({v3_slic3 = map(nc, supercell3, v3, 200, "euclidean", TRUE, "to_LAB")})
v3_slic3b = map2(v3_slic3, nc, add_cols)
v3_slic3c = do.call(rbind, v3_slic3b)

v3_rgb = st_rgb(stars::st_as_stars(raster::raster(v3))[,,,1:3],
                probs = c(0.02, 0.98), stretch = TRUE)

v3_p = ggplot() +
  geom_stars(data = v3_rgb)  +
  scale_fill_identity() +
  # scale_fill_continuous_sequential("Terrain 2") +
  geom_sf(data = v3_slic3c, fill = NA, color = "red", size = 0.4, aes(group = sizes)) +
  stars:::theme_stars() +
  transition_states(sizes, state_length = 2) +
  ease_aes("bounce-in") +
  enter_fade() +
  exit_shrink() +
  theme(legend.position = "none") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = NULL, y = NULL, title = "Compactness: {closest_state}")

# animate(g1a, renderer = ffmpeg_renderer())
anim_save("v3_p.gif", v3_p, height = 280/50, width = 550/50, units = "in", res = 150)

