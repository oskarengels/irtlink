# tests/testthat/test-ajk-haebara.R
# approximate jackknife for chain Haebara links, item estimating
# equations (Robitzsch, 2024, Stats, Eq. 14)

ajk_vs_jk_rel <- function(lk) {
  jk <- linking_error(lk, method = "jackknife")$le
  ajk <- linking_error(lk, method = "ajk")$le
  d <- c(ajk$le_mu[-1] - jk$le_mu[-1], ajk$le_sigma[-1] - jk$le_sigma[-1])
  ref <- pmax(abs(c(jk$le_mu[-1], jk$le_sigma[-1])), 1e-12)
  max(abs(d) / ref)
}

test_that("the Haebara estimating equations vanish at the sirt solution", {
  lk <- link_chain(ajk_hae_ipars(), method = "haebara")
  comp <- le_ipar_components(lk, unit_mode = "chain")
  st <- le_hae_settings(lk)
  Gi <- le_hae_gmat(le_hae_delta(lk), comp, st)
  G <- colSums(Gi)
  expect_lt(max(abs(G)) / mean(abs(Gi)), 1e-4)
})

test_that("the symmetric Haebara scores vanish at the sirt solution", {
  lk <- link_chain(ajk_hae_ipars(), method = "haebara", type = "symm")
  comp <- le_ipar_components(lk, unit_mode = "chain")
  st <- le_hae_settings(lk)
  Gi <- le_hae_gmat(le_hae_delta_variant(lk, st), comp, st)
  expect_lt(max(abs(colSums(Gi))) / mean(abs(Gi)), 1e-4)
})

test_that("the swapped-fit Haebara scores vanish for theta_metric = 'from'", {
  lk <- link_chain(ajk_hae_ipars(), method = "haebara",
                   theta_metric = "from")
  comp <- le_ipar_components(lk, unit_mode = "chain")
  st <- le_hae_settings(lk)
  Gi <- le_hae_gmat(le_hae_delta_variant(lk, st), comp, st)
  expect_lt(max(abs(colSums(Gi))) / mean(abs(Gi)), 1e-4)
})

test_that("AJK approximates the refit jackknife for chain Haebara", {
  for (case in list(list(pow = 2, weights = "uniform"),
                    list(pow = 2, weights = "normal_1"))) {
    lk <- link_chain(ajk_hae_ipars(), method = "haebara",
                     pow = case$pow, weights = case$weights)
    jk <- linking_error(lk, method = "jackknife")$le
    ajk <- linking_error(lk, method = "ajk")$le
    expect_identical(ajk$method, rep("ajk", nrow(ajk)))
    for (col in c("le_mu", "le_sigma")) {
      expect_equal(ajk[[col]][-1], jk[[col]][-1], tolerance = 0.02,
                   label = paste("haebara", case$weights, col))
    }
  }
})

test_that("the 3-group AJK converges to the refit jackknife as drift shrinks", {
  # One-step vs. fully refit delete-one: the relative deviation is of the
  # order of the drift (measured 6.5% / 1.7% / 0.7% at factors 1 / 0.25 /
  # 0.1), so it must both be small for moderate drift and shrink with it.
  rel_large <- ajk_vs_jk_rel(
    link_chain(ajk_hae_ipars(n_groups = 3), method = "haebara"))
  rel_small <- ajk_vs_jk_rel(
    link_chain(ajk_hae_ipars(n_groups = 3, drift = 0.25), method = "haebara"))
  expect_lt(rel_large, 0.10)
  expect_lt(rel_small, 0.03)
  expect_lt(rel_small, rel_large)
})

test_that("the analytic AJK is rejected internally for chain Haebara", {
  lk <- link_chain(ajk_hae_ipars(), method = "haebara")
  expect_error(le_ajk_chain(lk, method = "ajk", analytic = TRUE),
               "not available")
})

test_that("AJK for chain Haberman links is still rejected", {
  lk <- suppressWarnings(link_chain(ajk_hae_ipars(), method = "haberman"))
  expect_error(linking_error(lk, method = "ajk"), "currently available")
})
