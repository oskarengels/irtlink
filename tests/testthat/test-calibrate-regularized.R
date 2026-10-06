test_that("regularized (SBIC) recovers the trend on a warm-started eps grid", {
  skip_on_cran()
  set.seed(701)
  d <- sim_trend_data(n_groups = 3, N = 2500, I = 20, dif = "unbalanced",
                      dif_pct = 0.2, dif_effect = 0.7,
                      mu = c(0, 0.4, 0.8), sigma = c(1, 1, 1))
  cal <- calibrate(d, calibration = "regularized",
                   eps = c(0.01, 0.001, 0.0001))
  expect_identical(cal$calibration, "regularized")
  expect_true(cal$eps %in% c(0.01, 0.001, 0.0001))   # one eps was selected
  expect_equal(cal$trend$mu, c(0, 0.4, 0.8), tolerance = 0.15)
})

test_that("regularized (SBIC) accepts a single eps and uses it as-is", {
  skip_on_cran()
  set.seed(702)
  d <- sim_trend_data(n_groups = 3, N = 1800, I = 16, dif = "unbalanced",
                      dif_pct = 0.2, dif_effect = 0.7,
                      mu = c(0, 0.4, 0.8), sigma = c(1, 1, 1))
  cal <- calibrate(d, calibration = "regularized", eps = 0.01)
  expect_identical(cal$eps, 0.01)                    # no grid, the given eps
  expect_equal(cal$trend$mu, c(0, 0.4, 0.8), tolerance = 0.2)
})
