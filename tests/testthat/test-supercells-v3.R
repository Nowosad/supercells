v3_a = supercells(v3, 100, compactness = 1)
v3_b = supercells(v3, 100, compactness = 1, transform = "to_LAB")

test_that("supercells works for 3 var", {
  expect_equal(ncol(v3_a), 7)
  # expect_equal(nrow(v3_a), 80)
  # expect_equal(nrow(v3_b), 86)
  expect_equal(as.numeric(sf::st_bbox(v3_a)), unname(as.vector(terra::ext(v3)))[c(1, 3, 2, 4)])
})

# test_that("supercells matches reference (no geometry)", {
#   ref_path = testthat::test_path("testdata", "v3-supercells-v1.rds")
#   testthat::skip_if_not(file.exists(ref_path), "Reference file missing; create with old package version.")

#   current = list(
#     v3_a = sf::st_drop_geometry(v3_a),
#     v3_b = sf::st_drop_geometry(v3_b)
#   )
#   reference = readRDS(ref_path)
#   expect_equal(current, reference, tolerance = 1e-6)
# })
