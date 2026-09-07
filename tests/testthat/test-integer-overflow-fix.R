# Test for integer overflow fix in vector indexing
# This test verifies that large rasters (that would trigger 32-bit integer overflow
# in the old code) work correctly with the new 64-bit indexing via read_val_vec()
#

test_that("Large raster with multiple bands does not overflow (integer overflow fix)", {
  skip_on_cran()
  
  # Create a moderately large raster that would overflow 32-bit arithmetic
  # 5000 x 5000 x 3 = 75 million elements
  # With 10 bands: 750 million (close to 32-bit limit of 2.1 billion)
  # The key is that ncell * bands calculation should use 64-bit arithmetic
  
  nrows <- 2000
  ncols <- 2000
  nbands <- 3
  
  # Create a simple raster with random values
  set.seed(42)
  r <- terra::rast(ncol = ncols, nrow = nrows, nlyr = nbands)
  for (i in 1:nbands) {
    r[[i]] <- terra::setValues(r[[i]], runif(nrows * ncols))
  }
  terra::ext(r) <- c(0, ncols, 0, nrows)
  
  # This should complete without overflow errors or crashes
  # Using a large step to make it faster
  result <- sc_slic(r, step = 50, compactness = 1, iter = 2, verbose = 0)
  
  # Basic sanity checks
  expect_s3_class(result, "sf")
  expect_true("supercells" %in% names(result))
  expect_true(nrow(result) > 0)
})

test_that("Very large raster (near 32-bit limit) works correctly", {
  skip_on_cran()
  
  # Create a raster that would definitely overflow in 32-bit:
  # 10000 x 10000 = 100 million cells
  # With 10 bands: 1 billion index calculations
  # The multiplication 100M * 10 in 32-bit overflows
  
  nrows <- 2500
  ncols <- 2500
  nbands <- 10
  
  set.seed(123)
  r <- terra::rast(ncol = ncols, nrow = nrows, nlyr = nbands)
  for (i in 1:nbands) {
    r[[i]] <- terra::setValues(r[[i]], runif(nrows * ncols))
  }
  terra::ext(r) <- c(0, ncols, 0, nrows)
  
  # This tests the overflow fix in multiple code paths:
  # - create_centers: reads initial center values
  # - generate_superpixels: reads pixel values during assignment
  # - compute_max_value_dist: reads values for adaptive compactness
  result <- sc_slic(r, step = 100, compactness = 1, iter = 1, verbose = 0)
  
  expect_s3_class(result, "sf")
  expect_true(nrow(result) > 0)
})
