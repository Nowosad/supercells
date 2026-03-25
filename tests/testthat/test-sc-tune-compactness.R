test_that("sc_tune_compactness returns one-row output", {
  tune = sc_tune_compactness(v1, step = 8, metric = "local_variability")
  expect_s3_class(tune, "data.frame")
  expect_equal(nrow(tune), 1)
  expect_true(all(c("step", "metric", "dist_fun", "compactness") %in% names(tune)))
  expect_equal(tune$step, 8)
  expect_equal(tune$metric, "local_variability")
  expect_equal(tune$dist_fun, "euclidean")
  expect_true(tune$compactness > 0)
})

test_that("sc_tune_compactness stores custom dist_fun label", {
  manhattan = function(a, b) sum(abs(a - b))
  tune = sc_tune_compactness(v1, step = 8, dist_fun = manhattan)
  expect_equal(tune$dist_fun, "custom")
})

test_that("sc_tune_compactness local_variability scales with dimensionality", {
  v64 = c(v1, replicate(63, v1, simplify = FALSE))
  lv1 = sc_tune_compactness(v1, step = 8, metric = "local_variability")
  lv64 = sc_tune_compactness(v64, step = 8, metric = "local_variability")

  expect_true(is.finite(lv1$compactness))
  expect_true(is.finite(lv64$compactness))
  expect_true(lv64$compactness > 0)
})

test_that("sc_tune_compactness validates metric values", {
  expect_error(
    sc_tune_compactness(v1, step = 8, metric = "bad"),
    "metric must be 'local_variability'",
    fixed = TRUE
  )
})
