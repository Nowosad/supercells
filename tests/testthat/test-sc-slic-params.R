test_that("sc_slic_get_params/sc_slic_set_params work with minimal roundtrip", {
  sc = sc_slic(v1, step = 8, compactness = 1)
  params = sc_slic_get_params(sc)

  expect_s3_class(params, "data.frame")
  expect_equal(nrow(params), 1)
  expect_true(all(c("step", "compactness", "compactness_method", "dist_fun") %in% names(params)))
  expect_true(is.character(params$dist_fun))

  sc2 = sc
  attr(sc2, "step") = NULL
  attr(sc2, "compactness") = NULL
  attr(sc2, "compactness_method") = NULL
  attr(sc2, "dist_fun") = NULL
  sc2 = sc_slic_set_params(sc2, params)

  expect_equal(attr(sc2, "step"), attr(sc, "step"))
  expect_equal(attr(sc2, "compactness"), attr(sc, "compactness"))
  expect_equal(attr(sc2, "compactness_method"), attr(sc, "compactness_method"))
  expect_equal(attr(sc2, "dist_fun"), attr(sc, "dist_fun"))

  manhattan = function(a, b) sum(abs(a - b))
  sc_custom = sc_slic(v1, step = 8, compactness = 1, dist_fun = manhattan)
  params_custom = sc_slic_get_params(sc_custom)
  expect_true(is.na(params_custom$dist_fun))
})
