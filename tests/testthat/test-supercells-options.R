test_that("diagnostics are returned for single chunk", {
  v1_diag = supercells(v1, k = 30, compactness = 1, diagnostics = TRUE)
  diag = attr(v1_diag, "diagnostics")
  expect_true(is.list(diag))
  expect_true(all(c("iteration", "per_cluster", "per_pixel", "scale") %in% names(diag)))
})

# test_that("diagnostics are disabled when chunks are used", {
#   v1_diag_chunked = expect_warning(
#     supercells(v1, k = 30, compactness = 1, diagnostics = TRUE, chunks = 10),
#     "Diagnostics are only available"
#   )
#   expect_null(attr(v1_diag_chunked, "diagnostics"))
# })

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
