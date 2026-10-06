# tests/testthat/test-link-estimators.R

# Step truth for groups (t, t+1) with global mu/sigma (group 1 = N(0,1)):
# sigma_step = sigma[t+1] / sigma[t]; mu_step = (mu[t+1] - mu[t]) / sigma[t]

test_that("est_link_mm recovers step parameters exactly", {
  mu <- c(0, 0.3); sigma <- c(1, 1.1)
  ip <- as_ipars(make_identified_ipars(mu, sigma))
  est <- est_link_mm(ipars_pair(ip, 1, 2))
  expect_equal(est$sigma, 1.1, tolerance = 1e-8)
  expect_equal(est$mu, 0.3, tolerance = 1e-8)
})

test_that("est_link_mgm recovers step parameters exactly", {
  mu <- c(0, 0.3); sigma <- c(1, 1.1)
  ip <- as_ipars(make_identified_ipars(mu, sigma))
  est <- est_link_mgm(ipars_pair(ip, 1, 2))
  expect_equal(est$sigma, 1.1, tolerance = 1e-8)
  expect_equal(est$mu, 0.3, tolerance = 1e-8)
})

test_that("estimators handle a non-reference source group", {
  mu <- c(0, 0.3, 0.6); sigma <- c(1, 1.1, 1.2)
  ip <- as_ipars(make_identified_ipars(mu, sigma))
  est <- est_link_mgm(ipars_pair(ip, 2, 3))
  expect_equal(est$sigma, 1.2 / 1.1, tolerance = 1e-8)
  expect_equal(est$mu, (0.6 - 0.3) / 1.1, tolerance = 1e-8)
})

test_that("est_link_haberman recovers step truth", {
  mu <- c(0, 0.3); sigma <- c(1, 1.1)
  ip <- as_ipars(make_identified_ipars(mu, sigma))
  est <- est_link_haberman(ipars_pair(ip, 1, 2))
  expect_equal(est$sigma, 1.1, tolerance = 1e-4)
  expect_equal(est$mu, 0.3, tolerance = 1e-4)
})

test_that("estimators reject non-positive discriminations", {
  ip <- as_ipars(make_identified_ipars(c(0, 0.3), c(1, 1.1), I = 10))
  pars <- ipars_pair(ip, 1, 2)
  pars$a_from[1] <- 0
  expect_error(est_link_mm(pars), "Non-positive discrimination")
  expect_error(est_link_mgm(pars), "Non-positive discrimination")
  pars$a_from[1] <- -0.5
  expect_error(est_link_mm(pars), "Non-positive discrimination")
})

test_that("est_link_haebara and est_link_sl recover step truth", {
  mu <- c(0, 0.3); sigma <- c(1, 1.1)
  ip <- as_ipars(make_identified_ipars(mu, sigma))
  pars <- ipars_pair(ip, 1, 2)
  for (f in list(est_link_haebara, est_link_sl)) {
    est <- f(pars)
    expect_equal(est$sigma, 1.1, tolerance = 1e-4)
    expect_equal(est$mu, 0.3, tolerance = 1e-4)
  }
})

test_that("weight schemes do not change drift-free recovery", {
  ip <- as_ipars(make_identified_ipars(c(0, 0.3), c(1, 1.1)))
  pars <- ipars_pair(ip, 1, 2)
  for (w in c("uniform", "normal_0.5", "normal_1", "normal_2")) {
    est <- est_link_haebara(pars, weights = w)
    expect_equal(est$mu, 0.3, tolerance = 1e-4)
    expect_equal(est$sigma, 1.1, tolerance = 1e-4)
  }
})

test_that("pow = 0 path uses the robust eps rule without error", {
  ip <- as_ipars(make_identified_ipars(c(0, 0.3), c(1, 1.1)))
  est <- est_link_haebara(ipars_pair(ip, 1, 2), pow = 0)
  expect_equal(est$mu, 0.3, tolerance = 1e-3)
})

test_that("response function estimators reject non-positive discriminations and bad pow", {
  ip <- as_ipars(make_identified_ipars(c(0, 0.3), c(1, 1.1), I = 10))
  pars <- ipars_pair(ip, 1, 2)
  bad <- pars
  bad$a_from[1] <- -0.5
  expect_error(est_link_haebara(bad), "Non-positive discrimination")
  expect_error(est_link_sl(bad), "Non-positive discrimination")
  expect_error(est_link_haebara(pars, pow = 3), "pow")
})


test_that("est_link_haberman supports pow and use_intercepts", {
  ip <- as_ipars(make_identified_ipars(c(0, 0.3), c(1, 1.1)))
  pars <- ipars_pair(ip, 1, 2)
  for (p in c(2, 1, 0.5, 0)) {
    for (ui in c(TRUE, FALSE)) {
      est <- suppressWarnings(
        est_link_haberman(pars, pow = p, use_intercepts = ui)
      )
      expect_equal(est$mu, 0.3, tolerance = 1e-3,
                   label = sprintf("mu (pow=%s, use_intercepts=%s)", p, ui))
      expect_equal(est$sigma, 1.1, tolerance = 1e-3,
                   label = sprintf("sigma (pow=%s, use_intercepts=%s)", p, ui))
    }
  }
})

test_that("eps can be overridden without duplicate-argument errors", {
  ip <- as_ipars(make_identified_ipars(c(0, 0.3), c(1, 1.1)))
  pars <- ipars_pair(ip, 1, 2)
  est_h <- suppressWarnings(est_link_haberman(pars, eps = 0.01))
  expect_equal(est_h$mu, 0.3, tolerance = 1e-3)
  est_c <- est_link_haebara(pars, eps = 0.01)
  expect_equal(est_c$mu, 0.3, tolerance = 1e-3)
})

test_that("est_link_haberman validates pow like the response function estimators", {
  ip <- as_ipars(make_identified_ipars(c(0, 0.3), c(1, 1.1)))
  expect_error(est_link_haberman(ipars_pair(ip, 1, 2), pow = 3),
               "`pow` must be one of")
})

test_that("est_link_haberman rejects non-positive discriminations", {
  ip <- as_ipars(make_identified_ipars(c(0, 0.3), c(1, 1.1), I = 10))
  pars <- ipars_pair(ip, 1, 2)
  pars$a_from[1] <- -0.5
  expect_error(suppressWarnings(est_link_haberman(pars)),
               "Non-positive discrimination")
})
