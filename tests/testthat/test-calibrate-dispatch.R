test_that("calibrate rejects unknown engines and calibrations", {
  d <- sim_trend_data(n_groups = 2, N = 200, I = 14, dif = "none")
  expect_error(calibrate(d, engine = "tam"), "should be")
  expect_error(calibrate(d, engine = "xxirt"), "should be")
  expect_error(calibrate(d, calibration = "bogus"), "should be one of")
})

test_that("separate calibration leaves trend NULL; concurrent populates it", {
  skip_on_cran()
  set.seed(801)
  d <- sim_trend_data(n_groups = 2, N = 800, I = 14, dif = "none")
  expect_null(calibrate(d)$trend)
  cc <- calibrate(d, calibration = "concurrent")
  expect_false(is.null(cc$trend))
  expect_named(cc$trend, c("group", "mu", "sigma"))
})

test_that("print shows the trend when present", {
  skip_on_cran()
  set.seed(802)
  d <- sim_trend_data(n_groups = 2, N = 800, I = 14, dif = "none")
  cc <- calibrate(d, calibration = "concurrent")
  expect_output(print(cc), "Trend")
  # separate calibration must NOT print a Trend block
  sep <- calibrate(d)
  out <- capture.output(print(sep))
  expect_false(any(grepl("Trend", out)))
})

test_that("free_items is rejected for non-concurrent calibration", {
  d <- sim_trend_data(n_groups = 2, N = 200, I = 14, dif = "none")
  it <- colnames(as.data.frame(d[[1]]))[1]
  expect_error(calibrate(d, calibration = "separate", free_items = it),
               "concurrent")
})

test_that("data-based DIF diagnostics reject concurrent calibration", {
  skip_on_cran()
  set.seed(881)
  d <- sim_trend_data(n_groups = 2, N = 400, I = 12, dif = "none")
  cc <- calibrate(d, calibration = "concurrent")
  expect_error(responses_per_group(cc), "separate")
})
