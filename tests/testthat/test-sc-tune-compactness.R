test_that("sc_tune_compactness returns one-row data frame", {
  tune = sc_tune_compactness(v1, step = 8, iter = 1, sample_size = 500)
  expect_s3_class(tune, "data.frame")
  expect_equal(nrow(tune), 1)
  expect_true(all(c("step", "metric", "dist_fun", "compactness") %in% names(tune)))
  expect_equal(tune$step, 8)
  expect_equal(tune$metric, "global")
  expect_equal(tune$dist_fun, "euclidean")
  expect_true(tune$compactness > 0)
})

test_that("sc_tune_compactness returns local compactness when requested", {
  tune = sc_tune_compactness(v1, step = 8, iter = 1, sample_size = 500, metric = "local")
  expect_s3_class(tune, "data.frame")
  expect_equal(nrow(tune), 1)
  expect_true(all(c("step", "metric", "dist_fun", "compactness") %in% names(tune)))
  expect_equal(tune$step, 8)
  expect_equal(tune$metric, "local")
  expect_equal(tune$dist_fun, "euclidean")
  expect_true(tune$compactness > 0)
})

test_that("sc_tune_compactness stores custom dist_fun label", {
  manhattan = function(a, b) sum(abs(a - b))
  tune = sc_tune_compactness(v1, step = 8, iter = 1, sample_size = 500, dist_fun = manhattan)
  expect_equal(tune$dist_fun, "custom")
})
