test_that("sc_slic returns sf with attributes", {
  sc = sc_slic(v1, k = 30, compactness = 1)
  expect_s3_class(sc, "sf")
  expect_false(any(c("supercells", "x", "y") %in% names(sc)))
  expect_equal(attr(sc, "step"), attr(sc, "step"))
  expect_equal(attr(sc, "compactness"), 1)
  expect_equal(attr(sc, "method"), "slic")
  expect_true("supercells" %in% class(sc))
})

test_that("sc_slic handles step and centers", {
  set.seed(2021-11-21)
  centers = sf::st_sf(geom = sf::st_sample(sf::st_as_sfc(sf::st_bbox(v1)), 50, type = "random"))
  sc = sc_slic(v1, centers = centers, step = 8, compactness = 1, metadata = TRUE)
  expect_s3_class(sc, "sf")
})

test_that("sc_slic_raster returns raster with attributes", {
  sc_r = sc_slic_raster(v1, k = 30, compactness = 1)
  expect_s4_class(sc_r, "SpatRaster")
})

test_that("sc_slic_raster matches rasterized sc_slic", {
  sc = sc_slic(v1, k = 30, compactness = 1, metadata = TRUE)
  sc_r = sc_slic_raster(v1, k = 30, compactness = 1)
  ref_r = terra::rasterize(terra::vect(sc), v1, field = "supercells")
  ref_vals = terra::values(ref_r, mat = FALSE)
  out_vals = terra::values(sc_r, mat = FALSE)
  expect_true(isTRUE(all.equal(ref_vals, out_vals, check.attributes = FALSE)))
})

test_that("sc_slic validates arguments", {
  expect_error(sc_slic(v1, k = 10), "compactness", fixed = TRUE)
  expect_error(sc_slic(v1, k = 10, step = 5, compactness = 1), "either 'k' or 'step'", fixed = TRUE)
  expect_error(sc_slic(v1, centers = sf::st_sf(geom = sf::st_sfc()), compactness = 1), "step", fixed = TRUE)
})
