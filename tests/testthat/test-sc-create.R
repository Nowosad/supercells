test_that("sc_slic returns core output and attributes", {
  sc = sc_slic(v1, step = 8, compactness = 1)
  expect_s3_class(sc, "sf")
  expect_true("supercells" %in% class(sc))
  expect_true(all(c("supercells", "x", "y") %in% names(sc)))
  expect_false(is.null(attr(sc, "step")))
  expect_equal(attr(sc, "compactness"), 1)
  expect_equal(attr(sc, "compactness_method"), "constant")

  sc_auto = sc_slic(v1, step = 8, compactness = use_adaptive())
  expect_true(is.na(attr(sc_auto, "compactness")))
  expect_equal(attr(sc_auto, "compactness_method"), "local_max")
})

test_that("sc_slic supports custom centers", {
  set.seed(2021-11-21)
  centers = sf::st_sf(geom = sf::st_sample(sf::st_as_sfc(sf::st_bbox(v1)), 50, type = "random"))
  sc = sc_slic(v1, centers = centers, step = 8, compactness = 1)
  expect_s3_class(sc, "sf")
})

test_that("sc_slic does not retain assignments from earlier iterations", {
  set.seed(1)
  x = terra::rast(nrows = 5, ncols = 5, xmin = 0, xmax = 5,
                  ymin = 0, ymax = 5, crs = "EPSG:3857")
  terra::values(x) = runif(terra::ncell(x))
  centers = sf::st_as_sf(data.frame(x = runif(2, 0.5, 4.5),
                                    y = runif(2, 0.5, 4.5)),
                         coords = c("x", "y"), crs = 3857)

  iter1 = sc_slic_raster(x, centers = centers, step = 1,
                         compactness = 0.1, iter = 1, clean = FALSE)
  iter2 = sc_slic_raster(x, centers = centers, step = 1,
                         compactness = 0.1, iter = 2, clean = FALSE)

  expect_equal(terra::values(iter1, mat = FALSE)[[2]], 1)
  expect_equal(terra::values(iter2, mat = FALSE)[[2]], 2)
  expect_false(anyNA(terra::values(iter2, mat = FALSE)))
})

test_that("sc_slic_raster matches rasterized sc_slic output", {
  sc = sc_slic(v1, step = 8, compactness = 1, outcomes = "supercells")
  sc_r = sc_slic_raster(v1, step = 8, compactness = 1)
  expect_s4_class(sc_r, "SpatRaster")

  ref_r = terra::rasterize(terra::vect(sc), v1, field = "supercells")
  ref_vals = terra::values(ref_r, mat = FALSE)
  out_vals = terra::values(sc_r, mat = FALSE)
  expect_true(isTRUE(all.equal(ref_vals, out_vals, check.attributes = FALSE)))
})

test_that("sc_slic_raster assigns unique ids across chunks", {
  step = 8
  chunks = 10
  expect_warning(
    sc_r <- sc_slic_raster(v1, step = step, compactness = 1, chunks = chunks),
    "rounded up",
    fixed = TRUE
  )

  chunk_ext = .sc_chunk_extents(dim(v1), limit = ceiling(chunks / step) * step)
  ranges = lapply(seq_len(nrow(chunk_ext)), function(i) {
    ext = chunk_ext[i, ]
    chunk = sc_r[ext[1]:ext[2], ext[3]:ext[4], drop = FALSE]
    vals = terra::values(chunk, mat = FALSE)
    vals = vals[!is.na(vals)]
    if (length(vals) == 0) return(NULL)
    c(min = min(vals), max = max(vals))
  })
  ranges = Filter(Negate(is.null), ranges)

  if (length(ranges) > 1) {
    mins = vapply(ranges, function(x) x[["min"]], numeric(1))
    maxs = vapply(ranges, function(x) x[["max"]], numeric(1))
    expect_true(all(maxs[-length(maxs)] < mins[-1]))
  }
})

test_that("sc_slic assigns contiguous ids across chunks", {
  sc = sc_slic(v1, step = 8, compactness = 1, iter = 1, chunks = 16)

  expect_equal(sort(sc[["supercells"]]), seq_len(nrow(sc)))
})

test_that("sc_slic_raster enforces raster-specific guardrails", {
  expect_error(
    sc_slic_raster(v1, step = 8, compactness = 1, iter = 0),
    "iter = 0 returns centers only",
    fixed = TRUE
  )
  expect_error(
    sc_slic_raster(v1, step = 8, compactness = 1, outcomes = "values"),
    "supports only outcomes = 'supercells'",
    fixed = TRUE
  )
})

test_that("sc_slic_raster works with chunks = TRUE", {
  old_opt = getOption("supercells.chunk_mem_gb")
  options(supercells.chunk_mem_gb = 0.001)
  on.exit(options(supercells.chunk_mem_gb = old_opt), add = TRUE)

  sc_r = sc_slic_raster(v1, step = 8, compactness = 1, chunks = TRUE)
  expect_s4_class(sc_r, "SpatRaster")
  expect_equal(names(sc_r), "supercells")
})

test_that("sc_slic validates key argument errors", {
  expect_error(sc_slic(v1, k = 10), "compactness", fixed = TRUE)
  expect_error(sc_slic(v1, k = 0, compactness = 1), "single positive integer", fixed = TRUE)
  expect_error(sc_slic(v1, k = 1.5, compactness = 1), "single positive integer", fixed = TRUE)
  expect_error(sc_slic(v1, k = 10, step = 5, compactness = 1), "either 'k' or 'step'", fixed = TRUE)
  expect_error(sc_slic(v1, centers = sf::st_sf(geom = sf::st_sfc()), compactness = 1), "step", fixed = TRUE)
  expect_error(sc_slic(v1, step = 8, compactness = 1, dist_fun = "not_a_dist"), "does not exist", fixed = TRUE)
})

test_that("initial center values match their snapped raster cells", {
  x = terra::rast(nrows = 5, ncols = 5, xmin = 0, xmax = 5, ymin = 0, ymax = 5)
  terra::values(x) = seq_len(terra::ncell(x))

  pts = suppressWarnings(sc_slic_points(x, step = 4, compactness = 1, iter = 0))
  extracted = terra::extract(x, terra::vect(pts))[[2]]

  expect_equal(pts[[names(x)]], extracted)
})

test_that("custom centers use zero-based native cell indices", {
  x = terra::rast(nrows = 1, ncols = 1, xmin = 0, xmax = 1, ymin = 0, ymax = 1,
                  crs = "EPSG:3857")
  terra::values(x) = 5
  center = sf::st_as_sf(data.frame(x = 0.5, y = 0.5), coords = c("x", "y"), crs = 3857)

  pts = sc_slic_points(x, centers = center, step = 1, compactness = 1, iter = 0)

  expect_equal(nrow(pts), 1)
  expect_equal(pts[[names(x)]], 5)
})

test_that("sc_slic handles rasters with missing CRS", {
  x_no_crs = suppressWarnings({x = terra::rast(nrows = 40, ncols = 40, xmin = 0, xmax = 40, ymin = 0, ymax = 40); terra::crs(x) <- ""; terra::values(x) = seq_len(terra::ncell(x)); x})
  expect_no_warning(sc_slic(x_no_crs, k = 7, compactness = 30, avg_fun = median))
})

test_that("sc_slic handles step units and unit-related errors", {
  step_map = use_meters(8 * terra::res(v1)[1])
  sc = sc_slic(v1, step = step_map, compactness = 1)
  expect_s3_class(attr(sc, "step"), "units")
  expect_equal(as.numeric(units::drop_units(attr(sc, "step"))), as.numeric(units::drop_units(step_map)))

  x_lonlat = terra::rast(nrows = 50, ncols = 50, xmin = 0, xmax = 1, ymin = 0, ymax = 1, crs = "EPSG:4326")
  terra::values(x_lonlat) = seq_len(terra::ncell(x_lonlat))
  expect_error(
    sc_slic(x_lonlat, step = use_meters(1000), compactness = 1),
    "projected CRS",
    fixed = TRUE
  )

  expect_error(
    sc_slic(v1, step = units::set_units(1, "km"), compactness = 1),
    "use_meters(value)",
    fixed = TRUE
  )

  x_non_meter = terra::rast(nrows = 50, ncols = 50, xmin = 0, xmax = 5000, ymin = 0, ymax = 5000, crs = "EPSG:2277")
  terra::values(x_non_meter) = seq_len(terra::ncell(x_non_meter))
  expect_error(
    suppressWarnings(sc_slic(x_non_meter, step = use_meters(100), compactness = 1)),
    "meter units",
    fixed = TRUE
  )

  x_non_square = terra::rast(nrows = 10, ncols = 10, xmin = 0, xmax = 1000,
                             ymin = 0, ymax = 2000, crs = "EPSG:3857")
  terra::values(x_non_square) = seq_len(terra::ncell(x_non_square))
  expect_warning(
    sc_slic(x_non_square, step = use_meters(400), compactness = 1),
    "approximate for non-square cells",
    fixed = TRUE
  )
})

test_that("use_meters validates inputs", {
  x = use_meters(100)
  expect_s3_class(x, "units")
  expect_equal(as.character(units::deparse_unit(x)), "m")
  expect_error(use_meters(-1), "single positive number", fixed = TRUE)
  expect_error(use_meters(c(1, 2)), "single positive number", fixed = TRUE)
})

test_that("use_adaptive validates input", {
  x = use_adaptive()
  expect_s3_class(x, "sc_adaptive")
  expect_error(use_adaptive("other"), "must be 'local_max'", fixed = TRUE)
})
