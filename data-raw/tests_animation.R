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
# library(supercell)

# test 1 ------------------------------------------------------------------
srtm = rast(system.file("raster/srtm.tif", package = "spDataLarge"))
srtm_slic = supercell(srtm, 100, 1, "jensen_shannon")
ggplot() +
  geom_stars(data = stars::st_as_stars(raster::raster(srtm))) +
  scale_fill_continuous_sequential("Terrain 2") +
  geom_sf(data = srtm_slic, fill = NA, color = "red", size = 0.1) +
  stars:::theme_stars() +
  theme(legend.position = "none") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = NULL, y = NULL)


# jensen_shannon ----------------------------------------------------------
supercell2 = function(step, x, nc, dist_fun = "euclidean", clean = TRUE){
  supercell(x = x, step = step, nc = nc, dist_fun = dist_fun, clean = clean)
}
add_cols = function(x, col_value){
  x$sizes = as.factor(col_value)
  x
}

sizes = seq(20, 100, by = 20)
names(sizes) = sizes
system.time({srtm_slic1 = map(sizes, supercell2, srtm, 1, "jensen_shannon")})
srtm_slic1b = map2(srtm_slic1, sizes, add_cols)
srtm_slic1c = do.call(rbind, srtm_slic1b)

g1a = ggplot() +
  geom_stars(data = stars::st_as_stars(raster::raster(srtm))) +
  scale_fill_continuous_sequential("Terrain 2") +
  geom_sf(data = srtm_slic1c, fill = NA, color = "red", size = 0.1, aes(group = sizes)) +
  stars:::theme_stars() +
  transition_time(sizes) +
  ease_aes("linear") +
  enter_fade() +
  exit_shrink()

# animate(g1a, renderer = ffmpeg_renderer())
anim_save("g1a.gif", g1a)


# jensen_shannon2 ---------------------------------------------------------
supercell3 = function(nc, x, step, dist_fun = "euclidean", clean = TRUE){
  supercell(x = x, step = step, nc = nc, dist_fun = dist_fun, clean = clean)
}

nc = seq(1, 20, by = 2)
names(nc) = nc
system.time({srtm_slic2 = map(nc, supercell3, srtm, 20, "jensen_shannon")})
srtm_slic2b = map2(srtm_slic2, nc, add_cols)
srtm_slic2c = do.call(rbind, srtm_slic2b)

g2a = ggplot() +
  geom_stars(data = stars::st_as_stars(raster::raster(srtm))) +
  scale_fill_continuous_sequential("Terrain 2") +
  geom_sf(data = srtm_slic2c, fill = NA, color = "red", size = 0.1, aes(group = sizes)) +
  stars:::theme_stars() +
  transition_time(sizes) +
  ease_aes("linear") +
  enter_fade() +
  exit_shrink()

# animate(g1a, renderer = ffmpeg_renderer())
anim_save("g2a.gif", g2a)


# euclidean ---------------------------------------------------------------
nc = seq(400, 1, by = -15)
names(nc) = nc
system.time({srtm_slic3 = map(nc, supercell3, srtm, 20, "euclidean")})
srtm_slic3b = map2(srtm_slic3, nc, add_cols)
srtm_slic3c = do.call(rbind, srtm_slic3b)

g3a = ggplot() +
  geom_stars(data = stars::st_as_stars(raster::raster(srtm))) +
  scale_fill_continuous_sequential("Terrain 2") +
  geom_sf(data = srtm_slic3c, fill = NA, color = "#DA3DBE", size = 0.1, aes(group = sizes)) +
  stars:::theme_stars() +
  transition_states(sizes) +
  ease_aes("bounce-in") +
  enter_fade() +
  exit_shrink() +
  theme(legend.position = "none") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = NULL, y = NULL)

# animate(g1a, renderer = ffmpeg_renderer())
anim_save("g3a2.gif", g3a, height = 537/100, width = 435/100, units = "in", res = 150)

