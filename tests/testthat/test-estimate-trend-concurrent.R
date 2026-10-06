test_that("estimate_trend concurrent recovers the trend (joint/chain/restricted)", {
  skip_on_cran()
  set.seed(4201)
  d <- sim_trend_data(n_groups = 3, N = 2500, I = 18, dif = "none",
                      mu = c(0, 0.4, 0.8), sigma = c(1, 1.1, 1.2))
  for (ap in c("joint", "chain", "joint_restricted")) {
    tr <- estimate_trend(d, calibration = "concurrent", approach = ap)
    expect_equal(tr$trend$mu, c(0, 0.4, 0.8), tolerance = 0.12, info = ap)
    expect_equal(tr$trend$sigma, c(1, 1.1, 1.2), tolerance = 0.12, info = ap)
  }
})
