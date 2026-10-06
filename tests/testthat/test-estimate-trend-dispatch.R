test_that("fixed and regularized chain recover the trend", {
  skip_on_cran()
  set.seed(4401)
  d <- sim_trend_data(n_groups = 3, N = 2500, I = 20, dif = "none",
                      mu = c(0, 0.4, 0.8), sigma = c(1, 1.1, 1.2))
  trf <- estimate_trend(d, "fixed", "chain")
  trr <- estimate_trend(d, "regularized", "chain")
  expect_equal(trf$trend$mu, c(0, 0.4, 0.8), tolerance = 0.12)
  expect_equal(trr$trend$mu, c(0, 0.4, 0.8), tolerance = 0.15)
})

test_that("regularized joint estimation is rejected as unsupported", {
  set.seed(4403)
  d <- sim_trend_data(n_groups = 2, N = 300, I = 10, dif = "none")
  expect_error(estimate_trend(d, "regularized", "joint"), "not yet supported")
})

test_that("fixed ignores approach (FIPC is sequential; joint is not an error)", {
  skip_on_cran()
  set.seed(4402)
  d <- sim_trend_data(n_groups = 3, N = 1500, I = 14, dif = "none",
                      mu = c(0, 0.3, 0.6), sigma = c(1, 1, 1))
  trj <- estimate_trend(d, "fixed", "joint")     # no error; runs FIPC
  expect_s3_class(trj, "irtlink_trend")
  expect_equal(trj$trend$mu, c(0, 0.3, 0.6), tolerance = 0.12)
})

test_that("purify is rejected for fixed/regularized with a clear message", {
  d <- sim_trend_data(n_groups = 2, N = 300, I = 10, dif = "none")
  expect_error(estimate_trend(d, "regularized", "chain", dif = "purify"),
               "regulariz")
  expect_error(estimate_trend(d, "fixed", "chain", dif = "purify"),
               "not supported")
})
