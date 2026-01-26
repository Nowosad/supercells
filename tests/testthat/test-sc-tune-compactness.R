test_that("sc_tune_compactness returns one-row data frame", {
  tune = sc_tune_compactness(v1, step = 8, iter = 1, sample_size = 500)
  expect_s3_class(tune, "data.frame")
  expect_equal(nrow(tune), 1)
  expect_true(all(c("step", "median_value_dist", "median_spatial_dist", "compactness") %in% names(tune)))
  expect_equal(tune$step, 8)
  expect_true(tune$compactness > 0)
})
