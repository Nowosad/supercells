sc_full = sc_slic(v1, step = 8, compactness = 1,
                  outcomes = c("supercells", "coordinates", "values"))
sc_values = sc_slic(v1, step = 8, compactness = 1, outcomes = "values")
sc_auto = sc_slic(v1, step = 8, compactness = "auto",
                  outcomes = c("supercells", "coordinates", "values"))
manhattan = function(a, b) sum(abs(a - b))
sc_custom = sc_slic(v1, step = 8, compactness = 1, dist_fun = manhattan,
                    outcomes = c("supercells", "coordinates", "values"))

test_that("metrics outputs have expected structure", {
  pix = sc_metrics_pixels(v1, sc_full)
  expect_s4_class(pix, "SpatRaster")
  expect_equal(terra::nlyr(pix), 4)
  expect_true(all(c("spatial_scaled", "value_scaled", "combined", "balance") %in% names(pix)))

  cl = sc_metrics_supercells(v1, sc_full)
  expect_s3_class(cl, "sf")
  expect_true(all(c("supercells", "mean_value_dist_scaled", "mean_spatial_dist_scaled",
                    "mean_combined_dist", "balance") %in% names(cl)))

  gl = sc_metrics_global(v1, sc_full)
  expect_s3_class(gl, "data.frame")
  expect_equal(nrow(gl), 1)
  expect_true(all(c("step", "compactness", "n_supercells",
                    "mean_value_dist_scaled", "mean_spatial_dist_scaled",
                    "mean_combined_dist", "balance") %in% names(gl)))
})

test_that("metrics work without metadata columns", {
  expect_s4_class(sc_metrics_pixels(v1, sc_values), "SpatRaster")
  expect_s3_class(sc_metrics_supercells(v1, sc_values), "sf")
  expect_s3_class(sc_metrics_global(v1, sc_values), "data.frame")
})

test_that("metrics use stored attributes and dist_fun defaults", {
  gl = sc_metrics_global(v1, sc_full)
  expect_equal(gl$step, attr(sc_full, "step"))
  expect_equal(gl$compactness, attr(sc_full, "compactness"))

  gl_auto = sc_metrics_global(v1, sc_auto)
  expect_equal(gl_auto$compactness, "auto")

  g_attr = sc_metrics_global(v1, sc_custom)
  g_explicit = sc_metrics_global(v1, sc_custom, dist_fun = manhattan)
  expect_equal(g_attr$mean_value_dist_scaled, g_explicit$mean_value_dist_scaled, tolerance = 1e-8)
  expect_equal(g_attr$mean_spatial_dist_scaled, g_explicit$mean_spatial_dist_scaled, tolerance = 1e-8)
  expect_equal(g_attr$mean_combined_dist, g_explicit$mean_combined_dist, tolerance = 1e-8)
})

test_that("metrics follow step encoding in cells vs meters", {
  v1_map = terra::aggregate(v1, fact = 2, fun = mean, na.rm = TRUE)
  res_map = terra::res(v1_map)[1]
  step_cells = 8
  step_map = in_meters(step_cells * res_map)

  sc_cells = sc_slic(v1_map, step = step_cells, compactness = 1,
                     outcomes = c("supercells", "coordinates", "values"))
  sc_map = sc_slic(v1_map, step = step_map, compactness = 1,
                   outcomes = c("supercells", "coordinates", "values"))

  g_cells = sc_metrics_global(v1_map, sc_cells, scale = FALSE)
  g_map = sc_metrics_global(v1_map, sc_map, scale = FALSE)
  expect_equal(g_map$mean_spatial_dist / g_cells$mean_spatial_dist, res_map, tolerance = 1e-6)

  g_cells_scaled = sc_metrics_global(v1_map, sc_cells, scale = TRUE)
  g_map_scaled = sc_metrics_global(v1_map, sc_map, scale = TRUE)
  expect_equal(g_cells_scaled$mean_spatial_dist_scaled,
               g_map_scaled$mean_spatial_dist_scaled,
               tolerance = 1e-6)
})

test_that("metrics reject invalid dist_fun", {
  expect_error(sc_metrics_pixels(v1, sc_full, dist_fun = "not_a_dist"), "does not exist", fixed = TRUE)
})

test_that("metrics work after save/read when explicit args are supplied", {
  sc = sc_slic(v1, step = 8, compactness = 1,
               outcomes = c("supercells", "coordinates", "values"))
  gpkg = tempfile(fileext = ".gpkg")
  sf::st_write(sc, gpkg, quiet = TRUE)
  sc_disk = sf::st_read(gpkg, quiet = TRUE)
  unlink(gpkg)

  expect_error(
    sc_metrics_global(v1, sc_disk),
    "required when it is not stored in 'sc'",
    fixed = TRUE
  )

  gl = sc_metrics_global(v1, sc_disk, step = 8, compactness = 1, dist_fun = "euclidean")
  expect_s3_class(gl, "data.frame")
  expect_equal(nrow(gl), 1)

  pix = sc_metrics_pixels(
    v1, sc_disk,
    step = 8, compactness = 1, dist_fun = "euclidean",
    metrics = "combined"
  )
  expect_s4_class(pix, "SpatRaster")

  cl = sc_metrics_supercells(
    v1, sc_disk,
    step = 8, compactness = 1, dist_fun = "euclidean",
    metrics = "combined"
  )
  expect_s3_class(cl, "sf")
})
