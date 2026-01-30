test_that("sc_tune_compactness returns one-row data frame", {
  tune = sc_tune_compactness(v1, step = 8, iter = 1, sample_size = 500)
  expect_s3_class(tune, "data.frame")
  expect_equal(nrow(tune), 1)
  expect_true(all(c("step", "compactness_global") %in% names(tune)))
  expect_equal(tune$step, 8)
  expect_true(tune$compactness_global > 0)
})

test_that("sc_tune_compactness returns local compactness when requested", {
  tune = sc_tune_compactness(v1, step = 8, iter = 1, sample_size = 500, metrics = "local")
  expect_s3_class(tune, "data.frame")
  expect_equal(nrow(tune), 1)
  expect_true(all(c("step", "compactness_local") %in% names(tune)))
  expect_equal(tune$step, 8)
  expect_true(tune$compactness_local > 0)
})
