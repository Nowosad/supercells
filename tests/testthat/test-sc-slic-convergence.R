test_that("sc_slic_convergence returns plottable convergence data", {
  conv = sc_slic_convergence(v1, step = 8, compactness = 1, iter = 5)
  expect_s3_class(conv, "sc_slic_convergence")
  expect_s3_class(conv, "data.frame")
  expect_true(all(c("iter", "mean_distance") %in% names(conv)))
  expect_equal(nrow(conv), 5)
  expect_equal(conv$iter, seq_len(5))
  expect_true(all(is.finite(conv$mean_distance)))
  tmp = tempfile(fileext = ".pdf")
  grDevices::pdf(tmp)
  on.exit({
    grDevices::dev.off()
    unlink(tmp)
  }, add = TRUE)
  expect_silent(plot(conv))
})

test_that("sc_slic_convergence validates iter", {
  expect_error(
    sc_slic_convergence(v1, step = 8, compactness = 1, iter = 0),
    "iter must be > 0",
    fixed = TRUE
  )
})
