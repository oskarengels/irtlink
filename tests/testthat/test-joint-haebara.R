# tests/testthat/test-joint-haebara.R
# Package-internal joint Haebara linking: the simultaneous criterion with
# common item parameters (a_i, b_i) plus (mu_t, sigma_t), configurable
# quadrature weights, smoothed Lp loss. The reference-implementation
# comparison lives in the development oracle tests (tests_dev/).

test_that("joint Haebara recovers an exact drift-free trend", {
  mu <- c(0, 0.3, 0.6); sigma <- c(1, 1.1, 1.2)
  ip <- make_identified_ipars(mu, sigma, I = 10)
  for (w in c("normal_1", "uniform")) {
    lk <- link_joint(ip, method = "haebara", weights = w)
    expect_equal(lk$trend$mu, mu, tolerance = 1e-4, label = w)
    expect_equal(lk$trend$sigma, sigma, tolerance = 1e-4, label = w)
  }
})

test_that("the gradient vanishes at the fitted joint Haebara solution", {
  ip <- jh_ipars()
  for (w in c("uniform", "normal_0.5")) {
    lk <- link_joint(ip, method = "haebara", weights = w)
    fit <- lk$joint_fit
    g <- joint_hae_gradient(fit$par, fit$prep)
    expect_lt(max(abs(g)), 1e-4 * (1 + abs(joint_hae_value(fit$par,
                                                           fit$prep))))
  }
})

test_that("weight schemes change the drifted joint Haebara solution", {
  ip <- jh_ipars()
  lk_u <- link_joint(ip, method = "haebara", weights = "uniform")
  lk_n <- link_joint(ip, method = "haebara", weights = "normal_0.5")
  expect_true(all(is.finite(lk_u$trend$mu)))
  expect_true(all(is.finite(lk_n$trend$mu)))
  expect_gt(max(abs(lk_u$trend$mu - lk_n$trend$mu)), 1e-5)
})

test_that("the 1PL case fixes discriminations and scales", {
  mu <- c(0, 0.4)
  ip <- make_identified_ipars(mu, c(1, 1), I = 10)
  ip$a <- 1
  set.seed(5)
  ip$b[ip$group == 2] <- ip$b[ip$group == 2] + stats::rnorm(10, sd = 0.2)
  lk <- link_joint(ip, method = "haebara", weights = "uniform")
  expect_equal(lk$trend$sigma, c(1, 1), tolerance = 1e-10)
  expect_equal(lk$trend$mu[2], 0.4, tolerance = 0.1)
})

test_that("refit jackknife and restricted joint work with weights", {
  ip <- jh_ipars(n_groups = 2)
  lk <- link_joint(ip, method = "haebara", weights = "uniform")
  jk <- linking_error(lk, method = "jackknife")$le
  expect_gt(jk$le_mu[2], 0)
  lkr <- link_joint(jh_ipars(), method = "haebara", weights = "uniform",
                    restricted = TRUE)
  expect_identical(lkr$approach, "joint_restricted")
  expect_true(all(is.finite(lkr$trend$mu)))
})
