test_that("partial-invariance concurrent frees drifted items and recovers trend", {
  skip_on_cran()
  set.seed(501)
  d <- sim_trend_data(n_groups = 3, N = 2500, I = 20, dif = "unbalanced",
                      dif_pct = 0.2, dif_effect = 0.6,
                      mu = c(0, 0.4, 0.8), sigma = c(1, 1, 1))
  drift_items <- attr(d, "dif_items")
  expect_gt(length(drift_items), 0)
  cal <- calibrate(d, calibration = "concurrent", free_items = drift_items)
  expect_identical(cal$free_items, drift_items)
  expect_equal(cal$trend$mu, c(0, 0.4, 0.8), tolerance = 0.12)
})
