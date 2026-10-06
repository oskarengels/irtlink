# person sampling weights in the calibration strategies

sim_pw <- function(seed = 5, N = 400) {
  set.seed(seed)
  sim_trend_data(n_groups = 2, N = N, I = 10, overlap = 1, dif = "none")
}

test_that("unit weights reproduce the unweighted fit exactly", {
  d <- sim_pw()
  w1 <- lapply(d, function(x) rep(1, nrow(x)))
  for (cal_type in c("separate", "concurrent", "fixed")) {
    c0 <- calibrate(d, calibration = cal_type)
    cw <- calibrate(d, calibration = cal_type, pweights = w1)
    expect_identical(c0$ipars$a, cw$ipars$a)
    expect_identical(c0$ipars$b, cw$ipars$b)
  }
})

test_that("weight two equals duplicating every person", {
  d <- sim_pw(N = 250)
  w2 <- lapply(d, function(x) rep(2, nrow(x)))
  dd <- lapply(d, function(x) rbind(x, x))
  cw <- calibrate(d, pweights = w2, keep_vcov = TRUE)
  cd <- calibrate(dd, keep_vcov = TRUE)
  expect_equal(cw$ipars$a, cd$ipars$a, tolerance = 1e-6)
  expect_equal(cw$ipars$b, cd$ipars$b, tolerance = 1e-6)
  expect_equal(cw$vcov_ipars[[1]], cd$vcov_ipars[[1]], tolerance = 1e-4)
})

test_that("unequal weights move the estimates", {
  d <- sim_pw()
  wu <- lapply(d, function(x) rep(c(0.2, 1.8), length.out = nrow(x)))
  c0 <- calibrate(d)
  cw <- calibrate(d, pweights = wu)
  expect_gt(max(abs(c0$ipars$b - cw$ipars$b)), 1e-4)
})

test_that("weights work for GPCM and the regularized path", {
  set.seed(9)
  d <- sim_trend_data(n_groups = 2, N = 300, I = 8, overlap = 1,
                      model = "GPCM", n_cat = 3)
  w1 <- lapply(d, function(x) rep(1, nrow(x)))
  c0 <- calibrate(d, model = "GPCM")
  cw <- calibrate(d, model = "GPCM", pweights = w1)
  expect_identical(c0$ipars$b, cw$ipars$b)
  d2 <- sim_pw(N = 200)
  cr <- calibrate(d2, calibration = "regularized",
                  pweights = lapply(d2, function(x) rep(1, nrow(x))))
  expect_true(all(is.finite(cr$trend$mu)))
})

test_that("pweights validation and guards fire", {
  d <- sim_pw(N = 100)
  expect_error(calibrate(d, pweights = rep(1, 100)), "list")
  expect_error(calibrate(d, pweights = list(rep(1, 100))), "list")
  expect_error(calibrate(d, pweights = list(rep(1, 100), rep(-1, 100))),
               "positive")
  cw <- calibrate(d, pweights = lapply(d, function(x) rep(1, nrow(x))))
  expect_warning(detect_dif(cw), "person")
  expect_error(dependent_vcov(cw, ids = attr(d, "person_ids")),
               "person weights")
})

test_that("estimate_trend passes pweights through the dots", {
  d <- sim_pw()
  w1 <- lapply(d, function(x) rep(1, nrow(x)))
  t0 <- estimate_trend(d, link = "mgm")
  tw <- estimate_trend(d, link = "mgm", pweights = w1)
  expect_identical(t0$trend$mu, tw$trend$mu)
})
