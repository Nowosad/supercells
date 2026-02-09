test_that("sc_tune_compactness returns one-row output for both metrics", {
  for (metric in c("global", "local")) {
    tune = sc_tune_compactness(v1, step = 8, iter = 1, sample_size = 500, metric = metric)
    expect_s3_class(tune, "data.frame")
    expect_equal(nrow(tune), 1)
    expect_true(all(c("step", "metric", "dist_fun", "compactness") %in% names(tune)))
    expect_equal(tune$step, 8)
    expect_equal(tune$metric, metric)
    expect_equal(tune$dist_fun, "euclidean")
    expect_true(tune$compactness > 0)
  }
})

test_that("sc_tune_compactness stores custom dist_fun label", {
  manhattan = function(a, b) sum(abs(a - b))
  tune = sc_tune_compactness(v1, step = 8, iter = 1, sample_size = 500, dist_fun = manhattan)
  expect_equal(tune$dist_fun, "custom")
})

test_that("sc_tune_compactness value_scale rescales compactness and validates input", {
  g1 = sc_tune_compactness(
    v1, step = 8, iter = 1, metric = "global",
    value_scale = 1, sample_size = terra::ncell(v1)
  )
  g2 = sc_tune_compactness(
    v1, step = 8, iter = 1, metric = "global",
    value_scale = 2, sample_size = terra::ncell(v1)
  )
  expect_equal(g1$compactness / g2$compactness, 2, tolerance = 1e-6)

  l1 = sc_tune_compactness(v1, step = 8, iter = 1, metric = "local", value_scale = 1)
  l2 = sc_tune_compactness(v1, step = 8, iter = 1, metric = "local", value_scale = 2)
  expect_equal(l1$compactness / l2$compactness, 2, tolerance = 1e-6)

  expect_error(
    sc_tune_compactness(v1, step = 8, value_scale = 0),
    "value_scale must be a single positive number or 'auto'",
    fixed = TRUE
  )
  expect_error(
    sc_tune_compactness(v1, step = 8, value_scale = "bad"),
    "value_scale must be a single positive number or 'auto'",
    fixed = TRUE
  )
})
