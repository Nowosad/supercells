v1_a = supercells(v1, 100, compactness = 1)
v1_b = supercells(v1, 100, compactness = 1, clean = FALSE)
v1_c = supercells(v1, step = 8, compactness = 1)
my_minarea = 3
v1_d = supercells(v1, step = 8, compactness = 1, minarea = my_minarea)
v1_e = supercells(v1, step = 8, compactness = 0.1, avg_fun = "median", dist_fun = "jaccard")
v1_f = supercells(v1, 100, compactness = 1, avg_fun = mean)

test_that("supercells works for 1 var", {
  expect_equal(ncol(v1_a), 5)
  expect_equal(ncol(v1_a), ncol(v1_b))
  expect_equal(nrow(v1_a), 92)
  expect_equal(nrow(v1_b), 88)
  expect_equal(nrow(v1_e), 88)
  expect_equal(v1_a, v1_c)
  expect_true(all(as.numeric(sf::st_area(v1_d)) >= xres(v1) * yres(v1) * my_minarea))
  expect_equal(v1_a, v1_f)
})
