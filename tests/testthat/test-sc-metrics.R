test_that("sc_metrics_pixels returns raster with layers", {
  sc = sc_slic(v1, k = 30, compactness = 1, metadata = TRUE)
  pix = sc_metrics_pixels(sc, v1, compactness = 1, step = attr(sc, "step"))
  expect_s4_class(pix, "SpatRaster")
  expect_equal(terra::nlyr(pix), 3)
  expect_true(all(c("spatial", "value", "combined") %in% names(pix)))
})

test_that("sc_metrics_clusters returns sf with metrics", {
  sc = sc_slic(v1, k = 30, compactness = 1, metadata = TRUE)
  cl = sc_metrics_clusters(sc, v1, compactness = 1, step = attr(sc, "step"))
  expect_s3_class(cl, "sf")
  expect_true(all(c("supercells", "mean_value_dist", "mean_spatial_dist",
                    "mean_combined_dist", "compactness_ratio") %in% names(cl)))
})

test_that("sc_metrics_global returns single-row data.frame", {
  sc = sc_slic(v1, k = 30, compactness = 1, metadata = TRUE)
  gl = sc_metrics_global(sc, v1, compactness = 1, step = attr(sc, "step"))
  expect_s3_class(gl, "data.frame")
  expect_equal(nrow(gl), 1)
  expect_true(all(c("step", "compactness", "n_supercells",
                    "mean_value_dist", "mean_spatial_dist", "mean_combined_dist",
                    "compactness_ratio_mean", "mean_value_dist_w",
                    "mean_spatial_dist_w", "mean_combined_dist_w") %in% names(gl)))
})
