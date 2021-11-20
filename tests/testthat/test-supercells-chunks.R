v3_a = supercells(v3, 100, compactness = 1, chunk = 150)

test_that("supercells works for many chunks", {
  expect_equal(ncol(v3_a), 7)
  expect_equal(as.numeric(sf::st_bbox(v3_a)), unname(as.vector(terra::ext(v3)))[c(1, 3, 2, 4)])
})
