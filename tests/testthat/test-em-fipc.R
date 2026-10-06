test_that("fit_em_fipc recovers a trend and grows the bank", {
  skip_on_cran()
  set.seed(81)
  d <- sim_trend_data(n_groups = 3, N = 2000, I = 15, dif = "none",
                      mu = c(0, 0.3, 0.6), sigma = c(1, 1, 1),
                      design = "successive")
  fc <- fit_em_fipc(d, model = "2PL")
  expect_true(fc$converged)
  expect_equal(fc$trend$mu[1], 0)
  expect_equal(fc$trend$sigma[1], 1)
  expect_equal(fc$trend$mu, c(0, 0.3, 0.6), tolerance = 0.25)
  # every group's items are in the ipars table
  for (t in 1:3) {
    expect_setequal(fc$ipars$item[fc$ipars$group == t],
                    colnames(as.data.frame(d[[t]])))
  }
  expect_null(fc$model)
})

test_that("calibrate dispatches FIPC to the em engine", {
  skip_on_cran()
  set.seed(82)
  d <- sim_trend_data(n_groups = 2, N = 500, I = 15, dif = "none")
  cal <- calibrate(d, calibration = "fixed", engine = "em")
  expect_equal(cal$engine, "em")
  expect_equal(cal$trend$group, 1:2)
  expect_null(cal$models)
  expect_named(cal$converged, "fipc_fit")
})
