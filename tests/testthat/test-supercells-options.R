test_that("metadata columns can be removed", {
  v1_no_meta = supercells(v1, k = 30, compactness = 1, metadata = FALSE)
  expect_false(any(c("supercells", "x", "y") %in% names(v1_no_meta)))
})

test_that("invalid dist_fun errors", {
  expect_error(
    supercells(v1, k = 10, compactness = 1, dist_fun = "not_a_dist"),
    "does not exist",
    fixed = TRUE
  )
})
