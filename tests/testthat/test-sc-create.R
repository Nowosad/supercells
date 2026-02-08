test_that("sc_slic returns sf with attributes", {
  sc = sc_slic(v1, step = 8, compactness = 1)
  expect_s3_class(sc, "sf")
  expect_true(all(c("supercells", "x", "y") %in% names(sc)))
  expect_equal(attr(sc, "step"), attr(sc, "step"))
  expect_equal(attr(sc, "compactness"), 1)
  expect_true("supercells" %in% class(sc))
})

test_that("sc_slic stores adaptive compactness as 'auto'", {
  sc = sc_slic(v1, step = 8, compactness = "auto")
  expect_equal(attr(sc, "compactness"), "auto")
})

test_that("sc_slic handles step and centers", {
  set.seed(2021-11-21)
  centers = sf::st_sf(geom = sf::st_sample(sf::st_as_sfc(sf::st_bbox(v1)), 50, type = "random"))
  sc = sc_slic(v1, centers = centers, step = 8, compactness = 1,
               outcomes = c("supercells", "coordinates", "values"))
  expect_s3_class(sc, "sf")
})

test_that("sc_slic_raster returns raster with attributes", {
  sc_r = sc_slic_raster(v1, step = 8, compactness = 1)
  expect_s4_class(sc_r, "SpatRaster")
})

test_that("sc_slic_raster matches rasterized sc_slic", {
  sc = sc_slic(v1, step = 8, compactness = 1, outcomes = "supercells")
  sc_r = sc_slic_raster(v1, step = 8, compactness = 1)
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
    if (length(vals) == 0) {
      return(NULL)
    }
    c(min = min(vals), max = max(vals))
  })
  ranges = Filter(Negate(is.null), ranges)
  if (length(ranges) > 1) {
    mins = vapply(ranges, function(x) x[["min"]], numeric(1))
    maxs = vapply(ranges, function(x) x[["max"]], numeric(1))
    expect_true(all(maxs[-length(maxs)] < mins[-1]))
  }
})

test_that("auto chunk size aligns to step", {
  step = 8
  old_opt = getOption("supercells.chunk_mem_gb")
  options(supercells.chunk_mem_gb = 0.001)
  on.exit(options(supercells.chunk_mem_gb = old_opt), add = TRUE)
  wsize = .sc_chunk_optimize_size(dim(v1), getOption("supercells.chunk_mem_gb"), step = step)
  expect_true(wsize %% step == 0)
  expect_true(wsize >= step)
})

test_that("sc_slic validates arguments", {
  expect_error(sc_slic(v1, k = 10), "compactness", fixed = TRUE)
  expect_error(sc_slic(v1, k = 10, step = 5, compactness = 1), "either 'k' or 'step'", fixed = TRUE)
  expect_error(sc_slic(v1, centers = sf::st_sf(geom = sf::st_sfc()), compactness = 1), "step", fixed = TRUE)
  expect_error(sc_slic(v1, step = 8, compactness = 1, dist_fun = "not_a_dist"), "does not exist", fixed = TRUE)
})

test_that("sc_slic accepts units-based map step", {
  step_map = in_meters(8 * terra::res(v1)[1])
  sc = sc_slic(v1, step = step_map, compactness = 1)
  expect_s3_class(attr(sc, "step"), "units")
  expect_equal(as.numeric(units::drop_units(attr(sc, "step"))), as.numeric(units::drop_units(step_map)))
})

test_that("sc_slic rejects units-based step for lonlat rasters", {
  x = terra::rast(nrows = 50, ncols = 50, xmin = 0, xmax = 1, ymin = 0, ymax = 1, crs = "EPSG:4326")
  terra::values(x) = seq_len(terra::ncell(x))
  expect_error(
    suppressWarnings(sc_slic(x, step = in_meters(1000), compactness = 1)),
    "projected CRS",
    fixed = TRUE
  )
})

test_that("sc_slic rejects non-meter units for step", {
  expect_error(
    sc_slic(v1, step = units::set_units(1, "km"), compactness = 1),
    "must use meters",
    fixed = TRUE
  )
})

test_that("sc_slic rejects units-based step for projected non-meter CRS", {
  x = terra::rast(nrows = 50, ncols = 50, xmin = 0, xmax = 5000, ymin = 0, ymax = 5000, crs = "EPSG:2277")
  terra::values(x) = seq_len(terra::ncell(x))
  expect_error(
    suppressWarnings(sc_slic(x, step = in_meters(100), compactness = 1)),
    "meter units",
    fixed = TRUE
  )
})

test_that("in_meters validates inputs", {
  x = in_meters(100)
  expect_s3_class(x, "units")
  expect_equal(as.character(units::deparse_unit(x)), "m")
  expect_error(in_meters(-1), "single positive number", fixed = TRUE)
  expect_error(in_meters(c(1, 2)), "single positive number", fixed = TRUE)
})
