v1_a = supercells(v1, 100, compactness = 1)
v1_b = supercells(v1, 100, compactness = 1, clean = FALSE)

test_that("supercells works for 1 var", {
  expect_equal(ncol(v1_a), 5)
  expect_equal(ncol(v1_a), ncol(v1_b))
  expect_equal(nrow(v1_a), 112)
  expect_equal(nrow(v1_b), 108)
})
