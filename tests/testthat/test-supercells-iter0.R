set.seed(2021-11-20)
y = sf::st_sf(geom = sf::st_sample(sf::st_as_sfc(sf::st_bbox(v1)), 100, type = "random"))

vol_slic0a = supercells(v1, k = y, step = 10, compactness = 1, iter = 0)
vol_slic0b = supercells(v1, step = 10, compactness = 1, iter = 0)

test_that("supercells works for 0 iter", {
  expect_equal(ncol(vol_slic0a), 1)
  expect_equal(nrow(vol_slic0a), 100)
  expect_equal(nrow(vol_slic0b), 54)
})
