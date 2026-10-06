# tests/testthat/test-joint-pairwise.R
# Pairwise joint linking: one metric
# transformation per group, fitted by minimizing the summed IRF (Hae) or
# TCF (SL) discrepancy over ALL group pairs with common items. The
# reference-implementation comparison lives in the development oracle
# tests (tests_dev/).

test_that("pairwise joint recovers an exact drift-free trend", {
  mu <- c(0, 0.3, 0.6); sigma <- c(1, 1.1, 1.2)
  ip <- make_identified_ipars(mu, sigma, I = 10)
  for (m in c("haebara_pw", "sl_pw")) {
    lk <- link_joint(ip, method = m)
    expect_equal(lk$trend$mu, mu, tolerance = 1e-4, label = m)
    expect_equal(lk$trend$sigma, sigma, tolerance = 1e-4, label = m)
  }
})

test_that("pairwise joint works for a successive design without 1-3 overlap", {
  mu <- c(0, 0.3, 0.6); sigma <- c(1, 1.1, 1.2)
  ip <- make_identified_ipars(mu, sigma, I = 12)
  keep <- (ip$group == 1 & ip$item %in% sprintf("I%02d", 1:8)) |
          (ip$group == 2 & ip$item %in% sprintf("I%02d", 5:12)) |
          (ip$group == 3 & ip$item %in% sprintf("I%02d", 9:12))
  ip <- ip[keep, , drop = FALSE]
  lk <- link_joint(ip, method = "haebara_pw")
  expect_equal(lk$trend$mu, mu, tolerance = 1e-3)
  expect_equal(lk$trend$sigma, sigma, tolerance = 1e-3)
})

test_that("weight schemes change the drifted pairwise solution", {
  ip <- jp_ipars()
  lk_u <- link_joint(ip, method = "sl_pw", weights = "uniform")
  lk_n <- link_joint(ip, method = "sl_pw", weights = "normal_0.5")
  expect_true(all(is.finite(lk_u$trend$mu)))
  expect_gt(max(abs(lk_u$trend$mu - lk_n$trend$mu)), 1e-6)
})

test_that("the refit jackknife works for pairwise joint links", {
  ip <- jp_ipars(n_groups = 2)
  lk <- link_joint(ip, method = "haebara_pw", weights = "uniform")
  jk <- linking_error(lk, method = "jackknife")$le
  expect_gt(jk$le_mu[2], 0)
})
