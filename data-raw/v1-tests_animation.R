devtools::load_all()
library(spDataLarge)
library(terra)
library(sf)
library(stars)
library(ggspatial)
library(ggplot2)
library(gganimate)
library(gifski)
library(colorspace)
library(purrr)
# library(supercells)

# test 1 ------------------------------------------------------------------
srtm = rast(system.file("raster/srtm.tif", package = "spDataLarge"))
srtm_slic = supercells(srtm, 100, 1, "jsd")
ggplot() +
  geom_stars(data = stars::st_as_stars(srtm)) +
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
  geom_stars(data = stars::st_as_stars(srtm)) +
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
g3a_anim = animate(g3a, height = 537/100 * 150, width = 420/100 * 150, res = 150,
                   renderer = gifski_renderer())
anim_save("animation-compactness-srtm.gif", g3a_anim)


# rgb ---------------------------------------------------------------------
v3 = rast(system.file("raster/ortho.tif", package = "supercells"))

nc = c(seq(100, 50, by = -10), seq(50, 10, by = -2), seq(9, 1, by = -1), 0.5, 0.1)
names(nc) = nc
system.time({v3_slic3 = map(nc, supercell3, v3, 200, "euclidean", TRUE, "to_LAB")})
v3_slic3b = map2(v3_slic3, nc, add_cols)
v3_slic3c = do.call(rbind, v3_slic3b)

v3_rgb = st_rgb(stars::st_as_stars(v3)[,,,1:3],
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
v3_p_anim = animate(v3_p, height = 280/50 * 150, width = 550/50 * 150, res = 150,
                    renderer = gifski_renderer())
anim_save("animation-compactness-ortho.gif", v3_p_anim)


# tests iterations --------------------------------------------------------
supercell_i = function(iter, x, step, nc, dist_fun = "euclidean", clean = TRUE, transform = NULL){
  supercells(x = x, k = step, compactness = nc, dist_fun = dist_fun, clean = clean, iter = iter, transform = transform)
}
add_cols = function(x, col_value){
  x$sizes = as.factor(col_value)
  x
}
srtm = rast(system.file("raster/srtm.tif", package = "spDataLarge"))

i  = 1:10
names(i) = i
system.time({i_slic = map(i, supercell_i, srtm, 200, 50, clean = FALSE)})
i_slic2 = map2(i_slic, i, add_cols)
i_slic3 = do.call(rbind, i_slic2)

g_i = ggplot() +
  geom_stars(data = stars::st_as_stars(srtm)) +
  scale_fill_continuous_sequential("Terrain 2") +
  geom_sf(data = i_slic3, fill = NA, color = "#DA3DBE", size = 0.1, aes(group = sizes)) +
  stars:::theme_stars() +
  transition_states(sizes, state_length = 2) +
  ease_aes("bounce-in") +
  enter_fade() +
  exit_shrink() +
  theme(legend.position = "none") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = NULL, y = NULL, title = "Numer of interations: {closest_state}")

# animate(g1a, renderer = ffmpeg_renderer())
g_i_anim = animate(g_i, height = 537/100 * 150, width = 420/100 * 150, res = 150,
                   renderer = gifski_renderer())
anim_save("animation-iterations-srtm.gif", g_i_anim)

# tests iterations2 --------------------------------------------------------
v3 = rast(system.file("raster/ortho.tif", package = "supercells"))
i  = 1:30
names(i) = i
system.time({i_slic2 = map(i, supercell_i, v3, 200, 50, clean = FALSE, transform = "to_LAB")})
i_slic2 = map2(i_slic2, i, add_cols)
i_slic2 = do.call(rbind, i_slic2)

v3_rgb = st_rgb(stars::st_as_stars(v3)[,,,1:3],
                probs = c(0.02, 0.98), stretch = TRUE)

g_i2 = ggplot() +
  geom_stars(data = v3_rgb)  +
  scale_fill_identity() +
  # scale_fill_continuous_sequential("Terrain 2") +
  geom_sf(data = i_slic2, fill = NA, color = "red", size = 0.4, aes(group = sizes)) +
  stars:::theme_stars() +
  transition_states(sizes, state_length = 2) +
  ease_aes("bounce-in") +
  enter_fade() +
  exit_shrink() +
  theme(legend.position = "none") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = NULL, y = NULL, title = "Number of iterations: {closest_state}")

# animate(g1a, renderer = ffmpeg_renderer())
g_i2_anim = animate(g_i2, height = 280/50 * 150, width = 550/50 * 150, res = 150,
                    renderer = gifski_renderer())
anim_save("animation-iterations-ortho.gif", g_i2_anim)
