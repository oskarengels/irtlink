# tests/testthat/test-calibrate-integration.R
# Cross-strategy integration test: separate + chain linking, concurrent,
# and fixed/FIPC must all recover the same latent trend when there is no DIF.

test_that("all calibration strategies recover the same no-DIF trend", {
  skip_on_cran()
  set.seed(901)
  d <- sim_trend_data(n_groups = 3, N = 2500, I = 20, dif = "none",
                      mu = c(0, 0.4, 0.8), sigma = c(1, 1.1, 1.2))
  truth_mu <- c(0, 0.4, 0.8)
  sep  <- link_chain(calibrate(d), method = "mgm")$trend
  conc <- calibrate(d, calibration = "concurrent")$trend
  fipc <- calibrate(d, calibration = "fixed")$trend
  for (tr in list(sep, conc, fipc))
    expect_equal(tr$mu, truth_mu, tolerance = 0.12)
  truth_sigma <- c(1, 1.1, 1.2)
  for (tr in list(sep, conc, fipc))
    expect_equal(tr$sigma, truth_sigma, tolerance = 0.12)
})
