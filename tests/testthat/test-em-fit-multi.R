test_that("em_fit_multi with one group matches the separate fit", {
  skip_on_cran()
  set.seed(71)   # sim_trend_data needs >= 2 groups; simulate one directly
  J <- 12L; N <- 800L
  a <- runif(J, 0.6, 1.8); b <- rnorm(J)
  th <- rnorm(N)
  P <- plogis(sweep(outer(th, a), 2, a * b, "-"))
  X <- matrix(rbinom(N * J, 1L, P), N, J)
  colnames(X) <- paste0("i", seq_len(J))
  f1 <- em_fit_multi(list(X), "2PL")
  f0 <- em_fit(X, "2PL")
  m <- match(f0$item, f1$item)   # multi sorts the item union
  expect_equal(f1$a[m], f0$a, tolerance = 1e-6)
  expect_equal(f1$b[m], f0$b, tolerance = 1e-6)
  expect_equal(f1$trend, data.frame(group = 1L, mu = 0, sigma = 1))
})

test_that("em_fit_multi recovers a known trend under full invariance", {
  skip_on_cran()
  set.seed(72)
  d <- sim_trend_data(n_groups = 3, N = 3000, I = 15, dif = "none",
                      mu = c(0, 0.3, 0.6), sigma = c(1, 1.05, 1.1))
  fit <- em_fit_multi(d, "2PL")
  expect_true(fit$converged)
  expect_true(all(diff(fit$loglik_trace) > -1e-6))
  # seed 72's ML optimum sits at mu = (0, 0.245, 0.547): xxirt lands on
  # the same values, so the loose tolerance reflects sampling error, not
  # engine precision (the cross-engine gate covers precision).
  expect_equal(fit$trend$mu, c(0, 0.3, 0.6), tolerance = 0.2)
  expect_equal(fit$trend$sigma, c(1, 1.05, 1.1), tolerance = 0.2)
  expect_equal(fit$deviance, -2 * fit$loglik)
})

test_that("em_fit_multi partial invariance recovers a known drift", {
  skip_on_cran()
  set.seed(73)   # coherent simulation: one theta draw per person per group
  J <- 15L; N <- 6000L
  a <- runif(J, 0.6, 1.8); b <- rnorm(J)
  items <- paste0("i", seq_len(J))
  th1 <- rnorm(N, 0, 1)
  th2 <- rnorm(N, 0.3, 1)
  eta1 <- sweep(outer(th1, a), 2, a * b, "-")
  eta2 <- sweep(outer(th2, a), 2, a * b, "-")
  eta2[, 1] <- eta2[, 1] - 0.5           # uniform drift g = 0.5 on item 1
  X1 <- matrix(rbinom(N * J, 1L, plogis(eta1)), N, J)
  X2 <- matrix(rbinom(N * J, 1L, plogis(eta2)), N, J)
  colnames(X1) <- colnames(X2) <- items
  fit <- em_fit_multi(list(X1, X2), "2PL", free_items = "i1")
  expect_true("i1__G2_g" %in% names(fit$g_estimates))
  expect_equal(unname(fit$g_estimates["i1__G2_g"]), 0.5, tolerance = 0.2)
  expect_equal(fit$trend$mu[2], 0.3, tolerance = 0.2)
  expect_equal(fit$npar, 2L * J + 1L + 2L)
})

test_that("warm starts reproduce the cold-start optimum", {
  skip_on_cran()
  set.seed(76)
  d <- sim_trend_data(n_groups = 2, N = 1500, I = 12, dif = "none")
  cold <- em_fit_multi(d, "2PL")
  warm <- em_fit_multi(d, "2PL",
                       start = list(a = cold$a, nu = cold$nu, g = cold$g,
                                    mu = cold$trend$mu,
                                    sigma = cold$trend$sigma))
  expect_equal(warm$a, cold$a, tolerance = 1e-4)
  expect_equal(warm$b, cold$b, tolerance = 1e-4)
  expect_lte(warm$iter, cold$iter)
})

test_that("a penalized fit shrinks clean items and keeps drifted ones", {
  skip_on_cran()
  set.seed(77)   # coherent two-group simulation with one drifted item
  J <- 12L; N <- 3000L
  a <- runif(J, 0.7, 1.6); b <- rnorm(J)
  items <- paste0("i", seq_len(J))
  th1 <- rnorm(N, 0, 1); th2 <- rnorm(N, 0.3, 1)
  eta1 <- sweep(outer(th1, a), 2, a * b, "-")
  eta2 <- sweep(outer(th2, a), 2, a * b, "-")
  eta2[, 1] <- eta2[, 1] - 0.8
  X1 <- matrix(rbinom(N * J, 1L, plogis(eta1)), N, J)
  X2 <- matrix(rbinom(N * J, 1L, plogis(eta2)), N, J)
  colnames(X1) <- colnames(X2) <- items
  fit <- em_fit_multi(list(X1, X2), "2PL", free_items = items,
                      penalty = list(logN = log(2 * N), eps = 0.001))
  gs <- fit$g_estimates
  expect_gt(abs(gs["i1__G2_g"]), 0.4)
  expect_lt(max(abs(gs[setdiff(names(gs), "i1__G2_g")])), 0.15)
  expect_equal(fit$trend$mu[2], 0.3, tolerance = 0.15)
  expect_true(is.finite(fit$bic))
})

test_that("free_slope_items recovers nonuniform drift", {
  skip_on_cran()
  set.seed(78)
  J <- 12L; N <- 4000L
  a <- runif(J, 0.7, 1.6); b <- rnorm(J)
  items <- paste0("i", seq_len(J))
  th1 <- rnorm(N, 0, 1); th2 <- rnorm(N, 0.3, 1)
  eta1 <- sweep(outer(th1, a), 2, a * b, "-")
  eta2 <- sweep(outer(th2, a), 2, a * b, "-")
  eta2[, 1] <- eta2[, 1] + 0.5 * th2       # slope drift +0.5 on item 1
  X1 <- matrix(rbinom(N * J, 1L, plogis(eta1)), N, J)
  X2 <- matrix(rbinom(N * J, 1L, plogis(eta2)), N, J)
  colnames(X1) <- colnames(X2) <- items
  fit <- em_fit_multi(list(X1, X2), "2PL", free_slope_items = "i1")
  expect_true("i1__G2_h" %in% names(fit$h_estimates))
  expect_equal(unname(fit$h_estimates["i1__G2_h"]), 0.5, tolerance = 0.25)
  expect_length(fit$g_estimates, 0L)
  expect_equal(fit$npar, 2L * J + 1L + 2L)
})

test_that("the fit object carries the final expected counts", {
  skip_on_cran()
  set.seed(79)
  d <- sim_trend_data(n_groups = 2, N = 500, I = 12, dif = "none")
  fit <- em_fit_multi(d, "2PL")
  expect_length(fit$counts, 2L)
  for (t in 1:2) {
    expect_equal(dim(fit$counts[[t]]$njk),
                 c(length(fit$item), length(fit$theta)))
    expect_equal(sum(fit$counts[[t]]$nk),
                 nrow(as.data.frame(d[[t]])), tolerance = 1e-8)
  }
})

test_that("em_fit_multi 1PL fixes all discriminations", {
  skip_on_cran()
  set.seed(75)
  d <- sim_trend_data(n_groups = 2, N = 600, I = 12, dif = "none",
                      model = "1PL")
  fit <- em_fit_multi(d, "1PL")
  expect_true(all(fit$a == 1))
  expect_true(fit$converged)
})
