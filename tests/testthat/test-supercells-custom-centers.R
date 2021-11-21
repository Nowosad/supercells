set.seed(2021-11-21)
y1 = sf::st_sf(geom = sf::st_sample(sf::st_as_sfc(sf::st_bbox(v1)), 100, type = "random"))
y2 = sf::st_sf(geom = sf::st_sample(sf::st_as_sfc(sf::st_bbox(v1)), 100, type = "regular"))
y3 = sf::st_sf(geom = sf::st_sample(sf::st_as_sfc(sf::st_bbox(v1)), 100, type = "hexagonal"))

vol_slic1 = supercells(v1, k = y1, step = 10, compactness = 1, iter = 10, clean = TRUE)
vol_slic2 = supercells(v1, k = y2, step = 10, compactness = 1, iter = 10, clean = TRUE)
vol_slic3 = supercells(v1, k = y3, step = 10, compactness = 1, iter = 10, clean = TRUE)

test_that("supercells works for custom centers", {
  expect_equal(ncol(vol_slic1), 5)
  expect_equal(ncol(vol_slic2), 5)
  expect_equal(ncol(vol_slic3), 5)
})
