test_that("sc_metrics_pixels returns raster with layers", {
  sc = sc_slic(v1, step = 8, compactness = 1,
               outcomes = c("supercells", "coordinates", "values"))
  pix = sc_metrics_pixels(v1, sc)
  expect_s4_class(pix, "SpatRaster")
  expect_equal(terra::nlyr(pix), 4)
  expect_true(all(c("spatial_scaled", "value_scaled", "combined", "balance") %in% names(pix)))
})

test_that("sc_metrics functions work without metadata columns", {
  sc = sc_slic(v1, step = 8, compactness = 1, outcomes = "values")
  pix = sc_metrics_pixels(v1, sc)
  expect_s4_class(pix, "SpatRaster")

  cl = sc_metrics_supercells(v1, sc)
  expect_s3_class(cl, "sf")

  gl = sc_metrics_global(v1, sc)
  expect_s3_class(gl, "data.frame")
})

test_that("sc_metrics uses step and compactness from attributes", {
  sc = sc_slic(v1, step = 8, compactness = 1,
               outcomes = c("supercells", "coordinates", "values"))
  gl = sc_metrics_global(v1, sc)
  expect_equal(gl$step, attr(sc, "step"))
  expect_equal(gl$compactness, attr(sc, "compactness"))
})

test_that("sc_metrics uses adaptive compactness from compactness attribute", {
  sc = sc_slic(v1, step = 8, compactness = "auto",
               outcomes = c("supercells", "coordinates", "values"))
  gl = sc_metrics_global(v1, sc)
  expect_equal(gl$compactness, "auto")
})

test_that("sc_metrics_supercells returns sf with metrics", {
  sc = sc_slic(v1, step = 8, compactness = 1,
               outcomes = c("supercells", "coordinates", "values"))
  cl = sc_metrics_supercells(v1, sc)
  expect_s3_class(cl, "sf")
  expect_true(all(c("supercells", "mean_value_dist_scaled", "mean_spatial_dist_scaled",
                    "mean_combined_dist", "balance") %in% names(cl)))
  expect_true(all(is.finite(cl[["balance"]]) | is.na(cl[["balance"]])))
})

test_that("sc_metrics_global returns single-row data.frame", {
  sc = sc_slic(v1, step = 8, compactness = 1,
               outcomes = c("supercells", "coordinates", "values"))
  gl = sc_metrics_global(v1, sc)
  expect_s3_class(gl, "data.frame")
  expect_equal(nrow(gl), 1)
  expect_true(all(c("step", "compactness", "n_supercells",
                    "mean_value_dist_scaled", "mean_spatial_dist_scaled", "mean_combined_dist",
                    "balance") %in% names(gl)))
  expect_true(is.finite(gl[["balance"]]) | is.na(gl[["balance"]]))
})

test_that("sc_metrics invalid dist_fun errors", {
  sc = sc_slic(v1, step = 8, compactness = 1,
               outcomes = c("supercells", "coordinates", "values"))
  expect_error(sc_metrics_pixels(v1, sc, dist_fun = "not_a_dist"), "does not exist", fixed = TRUE)
})

test_that("sc_metrics spatial units follow step encoding", {
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

test_that("sc_metrics defaults to dist_fun attribute when missing", {
  manhattan = function(a, b) sum(abs(a - b))
  sc = sc_slic(v1, step = 8, compactness = 1, dist_fun = manhattan,
               outcomes = c("supercells", "coordinates", "values"))
  g_attr = sc_metrics_global(v1, sc)
  g_explicit = sc_metrics_global(v1, sc, dist_fun = manhattan)
  expect_equal(g_attr$mean_value_dist_scaled, g_explicit$mean_value_dist_scaled, tolerance = 1e-8)
  expect_equal(g_attr$mean_spatial_dist_scaled, g_explicit$mean_spatial_dist_scaled, tolerance = 1e-8)
  expect_equal(g_attr$mean_combined_dist, g_explicit$mean_combined_dist, tolerance = 1e-8)
})
