test_that("sc_metrics_pixels returns raster with layers", {
  sc = sc_slic(v1, step = 8, compactness = 1, metadata = TRUE)
  pix = sc_metrics_pixels(v1, sc)
  expect_s4_class(pix, "SpatRaster")
  expect_equal(terra::nlyr(pix), 3)
  expect_true(all(c("spatial_scaled", "value_scaled", "combined") %in% names(pix)))
})

test_that("sc_metrics functions work without metadata columns", {
  sc = sc_slic(v1, step = 8, compactness = 1, metadata = FALSE)
  pix = sc_metrics_pixels(v1, sc)
  expect_s4_class(pix, "SpatRaster")

  cl = sc_metrics_supercells(v1, sc)
  expect_s3_class(cl, "sf")

  gl = sc_metrics_global(v1, sc)
  expect_s3_class(gl, "data.frame")
})

test_that("sc_metrics uses step and compactness from attributes", {
  sc = sc_slic(v1, step = 8, compactness = 1, metadata = TRUE)
  gl = sc_metrics_global(v1, sc)
  expect_equal(gl$step, attr(sc, "step"))
  expect_equal(gl$compactness, attr(sc, "compactness"))
})

test_that("sc_metrics_supercells returns sf with metrics", {
  sc = sc_slic(v1, step = 8, compactness = 1, metadata = TRUE)
  cl = sc_metrics_supercells(v1, sc)
  expect_s3_class(cl, "sf")
  expect_true(all(c("supercells", "mean_value_dist_scaled", "mean_spatial_dist_scaled",
                    "mean_combined_dist", "balance") %in% names(cl)))
  expect_true(all(is.finite(cl[["balance"]]) | is.na(cl[["balance"]])))
})

test_that("sc_metrics_global returns single-row data.frame", {
  sc = sc_slic(v1, step = 8, compactness = 1, metadata = TRUE)
  gl = sc_metrics_global(v1, sc)
  expect_s3_class(gl, "data.frame")
  expect_equal(nrow(gl), 1)
  expect_true(all(c("step", "compactness", "n_supercells",
                    "mean_value_dist_scaled", "mean_spatial_dist_scaled", "mean_combined_dist",
                    "balance") %in% names(gl)))
  expect_true(is.finite(gl[["balance"]]) | is.na(gl[["balance"]]))
})

test_that("sc_metrics invalid dist_fun errors", {
  sc = sc_slic(v1, step = 8, compactness = 1, metadata = TRUE)
  expect_error(sc_metrics_pixels(v1, sc, dist_fun = "not_a_dist"), "does not exist", fixed = TRUE)
})
