test_that("fixed (FIPC) calibration recovers the true trend", {
  skip_on_cran()
  set.seed(601)
  d <- sim_trend_data(n_groups = 3, N = 2000, I = 20, dif = "none",
                      mu = c(0, 0.4, 0.8), sigma = c(1, 1.1, 1.2))
  cal <- calibrate(d, calibration = "fixed")
  expect_identical(cal$calibration, "fixed")
  expect_equal(cal$trend$mu,    c(0, 0.4, 0.8), tolerance = 0.12)
  expect_equal(cal$trend$sigma, c(1, 1.1, 1.2), tolerance = 0.12)
})

test_that("fixed 1PL calibration runs and fixes discrimination at 1", {
  skip_on_cran()
  set.seed(611)
  d <- sim_trend_data(n_groups = 3, N = 1500, I = 16, dif = "none",
                      model = "1PL", mu = c(0, 0.3, 0.6), sigma = c(1, 1, 1))
  cal <- calibrate(d, calibration = "fixed", model = "1PL")
  expect_true(all(abs(cal$ipars$a - 1) < 1e-6))
  expect_equal(cal$trend$mu, c(0, 0.3, 0.6), tolerance = 0.15)
})
