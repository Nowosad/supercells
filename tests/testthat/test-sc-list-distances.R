test_that("sc_list_distances returns the curated distance registry", {
  distances = sc_list_distances()
  expected = c(
    "euclidean", "squared_euclidean", "manhattan", "chebyshev", "canberra",
    "lorentzian", "soergel", "wavehedges", "czekanowski",
    "squared_chord", "squared_chi", "clark", "jsd", "non-intersection",
    "wasserstein1", "cosine_distance", "dtw", "dtw2d"
  )

  expect_s3_class(distances, "data.frame")
  expect_named(
    distances,
    c("name", "label", "input", "compares", "effect_of_layer_rescaling", "guidance")
  )
  expect_identical(distances$name, expected)
  expect_true(all(vapply(distances, is.character, logical(1))))
})

test_that("distance names are restricted to the curated set", {
  expect_error(
    .sc_util_prep_dist_fun("minkowski"),
    "'minkowski' is not supported",
    fixed = TRUE
  )
  expect_error(
    .sc_util_prep_dist_fun(c("euclidean", "manhattan")),
    "one non-empty string",
    fixed = TRUE
  )
})

test_that("every listed distance runs through sc_slic_raster", {
  x = terra::rast(nrows = 4, ncols = 4, nlyrs = 4, xmin = 0, xmax = 4,
                  ymin = 0, ymax = 4, crs = "EPSG:3857")
  values = seq_len(terra::ncell(x))
  input = cbind(values, rev(values), values %% 5 + 1, values %% 3 + 1)
  input = input / rowSums(input)
  terra::values(x) = input

  for (name in sc_list_distances()$name) {
    result = sc_slic_raster(x, step = 2, compactness = 1, dist_fun = name,
                            iter = 1, clean = FALSE)
    expect_false(anyNA(terra::values(result, mat = FALSE)), info = name)
  }
})

test_that("dtw2d requires paired coordinate layers", {
  x = terra::rast(nrows = 2, ncols = 2, nlyrs = 3, crs = "EPSG:3857")
  terra::values(x) = cbind(1:4, 2:5, 3:6)

  expect_error(
    sc_slic_raster(x, step = 2, compactness = 1, dist_fun = "dtw2d", clean = FALSE),
    "'dtw2d' requires an even number of layers",
    fixed = TRUE
  )
})

test_that("cosine_distance has defined zero-profile behavior", {
  x = terra::rast(nrows = 2, ncols = 2, nlyrs = 2, xmin = 0, xmax = 2,
                  ymin = 0, ymax = 2, crs = "EPSG:3857")
  terra::values(x) = rbind(c(0, 0), c(1e100, 2e100),
                           c(1e100, 2e100), c(1e100, 2e100))
  sc = sc_slic(x, step = 2, compactness = 1, dist_fun = "cosine_distance",
               iter = 1, clean = FALSE)
  value = sc_metrics_pixels(x, sc, metrics = "value", scale = FALSE)

  expect_equal(terra::values(value, mat = FALSE), c(1, 0, 0, 0), tolerance = 1e-12)
})

test_that("wasserstein1 accounts for distance between ordered bins", {
  x = terra::rast(nrows = 2, ncols = 2, nlyrs = 4, xmin = 0, xmax = 2,
                  ymin = 0, ymax = 2, crs = "EPSG:3857")
  terra::values(x) = rbind(c(1, 0, 0, 0), c(1, 0, 0, 0),
                           c(0, 0, 0, 1), c(0, 0, 0, 1))
  sc = sc_slic(x, step = 2, compactness = 1, dist_fun = "wasserstein1",
               iter = 1, clean = FALSE)
  value = sc_metrics_pixels(x, sc, metrics = "value", scale = FALSE)

  expect_equal(terra::values(value, mat = FALSE), rep(1.5, 4))
})
