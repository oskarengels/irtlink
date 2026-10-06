# tests/testthat/test-ajk-sl.R
# closed-form approximate jackknife for chain Stocking-Lord links
# (Robitzsch, 2025, Foundations 2, Eqs. 30-38)

sl_ajk_vs_jk_rel <- function(lk) {
  jk <- linking_error(lk, method = "jackknife")$le
  ajk <- linking_error(lk, method = "ajk")$le
  d <- c(ajk$le_mu[-1] - jk$le_mu[-1], ajk$le_sigma[-1] - jk$le_sigma[-1])
  ref <- pmax(abs(c(jk$le_mu[-1], jk$le_sigma[-1])), 1e-12)
  max(abs(d) / ref)
}

test_that("the SL AJK approximates the refit jackknife", {
  for (ws in c("uniform", "normal_2")) {
    lk <- link_chain(ajk_sl_ipars(), method = "sl", weights = ws)
    jk <- linking_error(lk, method = "jackknife")$le
    ajk <- linking_error(lk, method = "ajk")$le
    expect_identical(ajk$method, rep("ajk", nrow(ajk)))
    for (col in c("le_mu", "le_sigma")) {
      expect_equal(ajk[[col]][-1], jk[[col]][-1], tolerance = 0.05,
                   label = paste("sl", ws, col))
    }
  }
})

test_that("the SL AJK reproduces the reference implementation", {
  # Regression anchor: on this deterministic example the closed-form AJK
  # matched the article's reference implementation (lesl__0.156.R) to
  # all printed digits: le_mu = 0.07651, le_sigma = 0.03035 (2 groups,
  # uniform weights, verified 2026-09-20).
  lk <- link_chain(ajk_sl_ipars(), method = "sl")
  ajk <- linking_error(lk, method = "ajk")$le
  expect_equal(ajk$le_mu[2], 0.07651, tolerance = 2e-3)
  expect_equal(ajk$le_sigma[2], 0.03035, tolerance = 2e-3)
})

test_that("the 3-group SL AJK stays near the refit jackknife", {
  # Unlike the Haebara AJK, the paper's SL AJK does NOT converge to the
  # refit jackknife as the drift shrinks: Eq. 32 freezes D_I,t and drops
  # the D_t * dZ_i/d(delta) term, leaving a drift-independent relative
  # deviation (measured 3.7% to 6.7% across drift factors 1 to 0.1 at
  # I = 12; a full one-step with the complete delete-one score converges
  # to below 0.1%). The bound below documents that accuracy.
  rel_large <- sl_ajk_vs_jk_rel(
    link_chain(ajk_sl_ipars(n_groups = 3), method = "sl"))
  rel_small <- sl_ajk_vs_jk_rel(
    link_chain(ajk_sl_ipars(n_groups = 3, drift = 0.25), method = "sl"))
  expect_lt(rel_large, 0.10)
  expect_lt(rel_small, 0.10)
})

test_that("the chain SL AJK reports its closed form", {
  lk <- link_chain(ajk_sl_ipars(), method = "sl")
  le <- linking_error(lk, method = "ajk")$le
  expect_identical(le$method, rep("ajk", nrow(le)))
  expect_true(all(is.finite(le$le_mu)))
})

test_that("the SL AJK requires the quadratic loss", {
  lk <- link_chain(ajk_sl_ipars(), method = "sl", pow = 1)
  expect_error(linking_error(lk, method = "ajk"), "pow = 2")
})
