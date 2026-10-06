test_that("fit_em_sbic selects an eps and separates drifted items", {
  skip_on_cran()
  set.seed(83)
  J <- 12L; N <- 2500L
  a <- runif(J, 0.7, 1.6); b <- rnorm(J)
  items <- paste0("i", seq_len(J))
  th1 <- rnorm(N, 0, 1); th2 <- rnorm(N, 0.3, 1)
  eta1 <- sweep(outer(th1, a), 2, a * b, "-")
  eta2 <- sweep(outer(th2, a), 2, a * b, "-")
  eta2[, 2] <- eta2[, 2] - 0.8
  X1 <- matrix(rbinom(N * J, 1L, plogis(eta1)), N, J)
  X2 <- matrix(rbinom(N * J, 1L, plogis(eta2)), N, J)
  colnames(X1) <- colnames(X2) <- items
  fc <- fit_em_sbic(list(X1, X2), model = "2PL", eps = c(0.01, 0.001))
  expect_true(fc$eps %in% c(0.01, 0.001))
  expect_gt(abs(fc$g["i2__G2_g"]), 0.4)
  expect_lt(max(abs(fc$g[setdiff(names(fc$g), "i2__G2_g")])), 0.15)
  # seed 83's ML optimum sits at mu2 = 0.236 (xxirt SBIC: 0.235, and an
  # unpenalized oracle-partial fit gives the same value): the loose
  # tolerance reflects sampling error, the gate covers precision.
  expect_equal(fc$trend$mu[2], 0.3, tolerance = 0.25)
  expect_true(is.finite(fc$bic))
})

test_that("fit_em_sbic errors without common items", {
  skip_on_cran()
  X1 <- matrix(rbinom(40, 1L, 0.5), 20, 2,
               dimnames = list(NULL, c("a1", "a2")))
  X2 <- matrix(rbinom(40, 1L, 0.5), 20, 2,
               dimnames = list(NULL, c("b1", "b2")))
  expect_error(fit_em_sbic(list(X1, X2)), "no common items")
})

test_that("calibrate dispatches SBIC to the em engine", {
  skip_on_cran()
  set.seed(84)
  d <- sim_trend_data(n_groups = 2, N = 500, I = 15, dif = "none")
  cal <- calibrate(d, calibration = "regularized", engine = "em",
                   eps = 0.001)
  expect_equal(cal$engine, "em")
  expect_equal(cal$eps, 0.001)
  expect_equal(cal$trend$group, 1:2)
  expect_named(cal$converged, "sbic_fit")
})
