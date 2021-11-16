v3_a = supercells(v3, 100, compactness = 1)
v3_b = supercells(v3, 100, compactness = 1, transform = "to_LAB")

test_that("supercells works for 3 var", {
  expect_equal(ncol(v3_a), 7)
  expect_equal(nrow(v3_a), 79)
  expect_equal(nrow(v3_b), 80)
  expect_equal(as.numeric(st_bbox(v3_a)), as.numeric(st_bbox(v3)))
})
