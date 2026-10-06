# tests/testthat/test-sandwich-bc.R
# Bias-corrected sandwich linking errors with the item-diagonal meat
# correction B - D_tilde (Robitzsch, 2024, Stats 7, Eqs. 35/36). The SE
# component uses the full item-parameter covariance through
# U = -J_T A^{-1} H_gamma and never carries the finite-item factor.

sbc_methods <- c("sandwich_esw_bc", "sandwich_osw_bc", "sandwich_bosw_bc")

sbc_ipars <- function() {
  ip <- make_identified_ipars(c(0, 0.3, 0.6), c(1, 1.1, 1.2), I = 6)
  ip$b[ip$group == 2 & ip$item == "I01"] <-
    ip$b[ip$group == 2 & ip$item == "I01"] + 0.5
  ip$a[ip$group == 3 & ip$item == "I02"] <-
    ip$a[ip$group == 3 & ip$item == "I02"] * 1.2
  ip
}

sbc_psd_vcov <- function(lk, seed = 42) {
  n <- 2L * nrow(lk$ipars)
  set.seed(seed)
  M <- matrix(stats::rnorm(n * n, sd = 0.05), n, n)
  crossprod(M) / n
}

test_that("bc sandwich methods report the full decomposition and reduce LE", {
  lk <- link_chain(sbc_ipars(), method = "mgm")
  V <- sbc_psd_vcov(lk)
  for (method in sbc_methods) {
    le <- linking_error(lk, method = method, vcov = V)$le
    expect_named(le, c(
      "group", "le_mu", "le_sigma", "method", "n_reps",
      "se_mu", "se_sigma", "lebc_mu", "lebc_sigma",
      "te_mu", "te_sigma", "tebc_mu", "tebc_sigma"
    ))
    expect_identical(le$method, rep(method, nrow(le)))
    expect_true(all(is.finite(le$lebc_mu)))
    expect_gt(max(le$se_mu), 0)
    # D_tilde is positive semidefinite, so the correction cannot increase LE.
    expect_true(all(le$lebc_mu <= le$le_mu + 1e-12))
    expect_true(all(le$lebc_sigma <= le$le_sigma + 1e-12))
  }
})

test_that("with a zero vcov the bc sandwich equals the plain sandwich", {
  lk <- link_chain(sbc_ipars(), method = "mgm")
  V0 <- matrix(0, nrow = 2L * nrow(lk$ipars), ncol = 2L * nrow(lk$ipars))
  plain_of <- c(sandwich_esw_bc = "sandwich_esw",
                sandwich_osw_bc = "sandwich_osw",
                sandwich_bosw_bc = "sandwich_bosw")
  for (method in sbc_methods) {
    bc <- linking_error(lk, method = method, vcov = V0)$le
    plain <- linking_error(lk, method = plain_of[[method]])$le
    expect_equal(bc$se_mu, rep(0, nrow(bc)))
    expect_equal(bc$se_sigma, rep(0, nrow(bc)))
    expect_equal(bc$lebc_mu, bc$le_mu, tolerance = 1e-12)
    expect_equal(bc$lebc_sigma, bc$le_sigma, tolerance = 1e-12)
    # Analytic bread/derivatives vs. the numeric ones of the plain methods.
    expect_equal(bc$le_mu, plain$le_mu, tolerance = 1e-6)
    expect_equal(bc$le_sigma, plain$le_sigma, tolerance = 1e-6)
  }
})

test_that("bosw_bc carries the finite-item factor on LE but not on SE", {
  lk <- link_chain(sbc_ipars(), method = "mgm")
  V <- sbc_psd_vcov(lk)
  osw <- linking_error(lk, method = "sandwich_osw_bc", vcov = V)$le
  bosw <- linking_error(lk, method = "sandwich_bosw_bc", vcov = V)$le
  n_units <- unique(osw$n_reps)
  fac <- sqrt(n_units / (n_units - 1))
  expect_equal(bosw$le_mu, fac * osw$le_mu, tolerance = 1e-12)
  expect_equal(bosw$le_sigma, fac * osw$le_sigma, tolerance = 1e-12)
  expect_equal(bosw$lebc_mu, fac * osw$lebc_mu, tolerance = 1e-12)
  expect_equal(bosw$lebc_sigma, fac * osw$lebc_sigma, tolerance = 1e-12)
  expect_equal(bosw$se_mu, osw$se_mu, tolerance = 1e-12)
  expect_equal(bosw$se_sigma, osw$se_sigma, tolerance = 1e-12)
})

test_that("bc sandwich SE equals the analytic jackknife_bc SE", {
  lk <- link_chain(sbc_ipars(), method = "mgm")
  V <- sbc_psd_vcov(lk)
  sw <- linking_error(lk, method = "sandwich_osw_bc", vcov = V)$le
  jk <- linking_error(lk, method = "jackknife_bc", vcov = V)$le
  expect_equal(sw$se_mu, jk$se_mu, tolerance = 1e-10)
  expect_equal(sw$se_sigma, jk$se_sigma, tolerance = 1e-10)
})

test_that("bc sandwich methods work for joint PHL links", {
  for (iw in c("uniform", "inverse_admin")) {
    lk <- link_joint(sbc_ipars(), method = "phl", item_weights = iw)
    V <- sbc_psd_vcov(lk)
    le <- linking_error(lk, method = "sandwich_osw_bc", vcov = V)$le
    expect_gt(max(le$le_mu), 0)
    expect_true(all(is.finite(le$lebc_mu)))
    expect_true(all(le$lebc_mu <= le$le_mu + 1e-12))
    plain <- linking_error(lk, method = "sandwich_osw")$le
    expect_equal(le$le_mu, plain$le_mu, tolerance = 1e-6)
    expect_equal(le$le_sigma, plain$le_sigma, tolerance = 1e-6)
  }
})

test_that("bc sandwich methods require vcov", {
  lk <- link_chain(sbc_ipars(), method = "mgm")
  for (method in sbc_methods) {
    expect_error(linking_error(lk, method = method), "vcov")
  }
})

# Regression tests for the dependent-samples D_tilde: the correction keeps
# the full within-item blocks across groups (cross-group covariances of the
# same item included), not only the per-(group, item) 2x2 blocks.

test_that("D_tilde keeps within-item cross-group covariance blocks", {
  lk <- link_chain(sbc_ipars(), method = "mgm")
  ip <- lk$ipars
  V <- sbc_psd_vcov(lk)
  # item-block-diagonal version of V (within-item blocks across groups),
  # rows in the order of lk$ipars, two parameters per (group, item) row
  W <- V * 0
  for (it in unique(ip$item)) {
    rows <- which(ip$item == it)
    idx <- as.integer(rbind(2L * rows - 1L, 2L * rows))
    W[idx, idx] <- V[idx, idx]
  }
  # 2x2 group-item block-diagonal version (the old, incorrect restriction)
  W2 <- V * 0
  for (r in seq_len(nrow(ip))) {
    idx <- c(2L * r - 1L, 2L * r)
    W2[idx, idx] <- V[idx, idx]
  }
  le_full <- linking_error(lk, method = "sandwich_osw_bc", vcov = V)$le
  le_wblk <- linking_error(lk, method = "sandwich_osw_bc", vcov = W)$le
  le_gblk <- linking_error(lk, method = "sandwich_osw_bc", vcov = W2)$le
  # The subtracted projection is the SE built from the item-block-diagonal
  # covariance: lebc^2 = le^2 - se(W)^2 for every group and parameter.
  for (col in c("mu", "sigma")) {
    lebc2 <- pmax(le_full[[paste0("le_", col)]]^2 -
                    le_wblk[[paste0("se_", col)]]^2, 0)
    expect_equal(le_full[[paste0("lebc_", col)]], sqrt(lebc2),
                 tolerance = 1e-10)
  }
  # The correction must depend on the within-item cross-group blocks:
  # W and W2 differ only in those blocks and must give different lebc.
  expect_false(isTRUE(all.equal(le_wblk$lebc_mu, le_gblk$lebc_mu,
                                tolerance = 1e-8)))
  # Nothing outside the within-item blocks is used for D_tilde: the full
  # covariance and its item-block-diagonal restriction give the same lebc.
  expect_equal(le_full$lebc_mu, le_wblk$lebc_mu, tolerance = 1e-12)
  expect_equal(le_full$lebc_sigma, le_wblk$lebc_sigma, tolerance = 1e-12)
})

test_that("te and tebc are consistent with the truncated le columns", {
  lk <- link_chain(sbc_ipars(), method = "mgm")
  # inflate the covariance so that the corrected LE variance is negative
  V <- sbc_psd_vcov(lk) * 9
  for (method in c(sbc_methods, "jackknife_bc")) {
    le <- linking_error(lk, method = method, vcov = V)$le
    expect_equal(le$te_mu, sqrt(le$se_mu^2 + le$le_mu^2), tolerance = 1e-12)
    expect_equal(le$te_sigma, sqrt(le$se_sigma^2 + le$le_sigma^2),
                 tolerance = 1e-12)
    expect_equal(le$tebc_mu, sqrt(le$se_mu^2 + le$lebc_mu^2),
                 tolerance = 1e-12)
    expect_equal(le$tebc_sigma, sqrt(le$se_sigma^2 + le$lebc_sigma^2),
                 tolerance = 1e-12)
    # truncation happens on the LE component, never on the sum: the
    # bias-corrected TE is at least the SE.
    expect_true(all(le$tebc_mu >= le$se_mu - 1e-12))
    expect_true(all(le$tebc_sigma >= le$se_sigma - 1e-12))
  }
})
