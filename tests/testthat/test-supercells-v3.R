v3_a = supercells(v3, 100, compactness = 1)
v3_b = supercells(v3, 100, compactness = 1, transform = "to_LAB")

test_that("supercells works for 3 var", {
  expect_equal(ncol(v3_a), 7)
  expect_equal(nrow(v3_a), 83)
  expect_equal(nrow(v3_b), 93)
  expect_equal(as.numeric(sf::st_bbox(v3_a)), unname(as.vector(terra::ext(v3)))[c(1, 3, 2, 4)])
})

