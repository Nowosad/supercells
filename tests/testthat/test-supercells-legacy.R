test_that("legacy supercells wrapper handles basic options", {
  sc = supercells(v1, k = 100, compactness = 1)
  expect_s3_class(sc, "sf")
  expect_true(all(c("supercells", "x", "y") %in% names(sc)))

  sc_no_meta = supercells(v1, k = 30, compactness = 1, metadata = FALSE)
  expect_false(any(c("supercells", "x", "y") %in% names(sc_no_meta)))

  expect_error(
    supercells(v1, k = 10, compactness = 1, dist_fun = "not_a_dist"),
    "does not exist",
    fixed = TRUE
  )
})

test_that("legacy supercells wrapper supports iter = 0 and custom centers", {
  set.seed(2021-11-20)
  centers = sf::st_sf(geom = sf::st_sample(sf::st_as_sfc(sf::st_bbox(v1)), 100, type = "random"))

  sc_custom = supercells(v1, k = centers, step = 10, compactness = 1, iter = 0)
  sc_default = supercells(v1, step = 10, compactness = 1, iter = 0)
  expect_s3_class(sc_custom, "sf")
  expect_s3_class(sc_default, "sf")
  expect_true(nrow(sc_custom) > 0)
  expect_true(nrow(sc_default) > 0)
})

test_that("legacy supercells wrapper supports transform and chunked calls", {
  sc_lab = supercells(v3, 100, compactness = 1, transform = "to_LAB")
  sc_chunk = supercells(v3, 100, compactness = 1, chunk = 150)

  expect_s3_class(sc_lab, "sf")
  expect_s3_class(sc_chunk, "sf")
  expect_equal(as.numeric(sf::st_bbox(sc_lab)), unname(as.vector(terra::ext(v3)))[c(1, 3, 2, 4)])
  expect_equal(as.numeric(sf::st_bbox(sc_chunk)), unname(as.vector(terra::ext(v3)))[c(1, 3, 2, 4)])
})
