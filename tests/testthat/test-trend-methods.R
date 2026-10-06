make_trend_fixture <- function() {
  structure(list(
    trend = data.frame(group = 1:3, mu = c(0, 0.4, 0.8), sigma = c(1, 1.1, 1.2),
                       se = NA_real_, le = NA_real_, te = NA_real_),
    calibration = "separate", approach = "chain", link = "mgm",
    calib = NULL, dif = NULL, link_obj = NULL, call = NULL),
    class = "irtlink_trend")
}
test_that("print/summary/coef work for irtlink_trend", {
  tr <- make_trend_fixture()
  expect_output(print(tr), "irtlink trend")
  expect_output(print(tr), "separate / chain")
  expect_identical(coef(tr), tr$trend)
})
test_that("plot draws a trend without error", {
  tr <- make_trend_fixture()
  grDevices::pdf(NULL); on.exit(grDevices::dev.off())
  expect_no_error(plot(tr))
})

test_that("print/summary handle an irtlink_calib with NULL ipars", {
  # composed concurrent fits (chain / restricted) embed a calib without a
  # single invariant ipars set; the S3 methods must not crash on it.
  cal <- structure(
    list(ipars = NULL, model = "2PL", calibration = "concurrent",
         engine = NA_character_, n_groups = 3L, N = c(100L, 100L, 100L),
         trend = data.frame(group = 1:3, mu = c(0, 0.4, 0.8),
                            sigma = c(1, 1, 1)), converged = TRUE),
    class = "irtlink_calib")
  expect_no_error(print(cal))
  expect_no_error(summary(cal))
  expect_null(coef(cal))
})
