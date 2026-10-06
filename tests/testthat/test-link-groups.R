# link_groups(), selectable reference group, and the new simulator
# options (non-uniform DIF, dependent person samples).

toy_ipars <- function() {
  mu_t <- c(0, 0.3, 0.6); sigma_t <- c(1, 1.1, 1.2)
  b_ref <- seq(-1.5, 1.5, length.out = 10)
  data.frame(
    group = rep(1:3, each = 10),
    item = rep(sprintf("I%02d", 1:10), 3),
    a = rep(sigma_t, each = 10),
    b = c(b_ref, (b_ref - 0.3) / 1.1, (b_ref - 0.6) / 1.2)
  )
}

test_that("link_groups dispatches to link_chain and link_joint", {
  ip <- toy_ipars()
  lg <- link_groups(ip, approach = "chain", method = "mgm")
  lc <- link_chain(ip, method = "mgm")
  expect_equal(lg$trend, lc$trend)
  expect_identical(lg$approach, "chain")

  lgj <- link_groups(ip, approach = "joint", method = "haberman",
                     variant = "pairwise")
  lj <- link_joint(ip, method = "phl")
  expect_equal(lgj$trend, lj$trend)
  expect_identical(lgj$approach, "joint")

  lgr <- link_groups(ip, approach = "joint_restricted",
                     method = "haberman", variant = "pairwise")
  expect_identical(lgr$approach, "joint_restricted")
})

test_that("link_groups rejects approach-incompatible choices", {
  ip <- toy_ipars()
  expect_error(link_groups(ip, approach = "chain", method = "phl"),
               "must be one of")
  expect_error(link_groups(ip, approach = "joint", method = "mgm"),
               "must be one of")
  expect_error(link_groups(ip, approach = "chain", method = "mgm",
                           variant = "pairwise"),
               "pairwise")
})

test_that("variant = pairwise maps onto the pairwise back ends", {
  ip <- toy_ipars()
  lp <- link_groups(ip, approach = "joint", method = "haberman",
                    variant = "pairwise")
  expect_identical(lp$method, "phl")
  expect_equal(lp$trend, link_joint(ip, method = "phl")$trend)
  lh <- link_groups(ip, approach = "joint", method = "haebara",
                    variant = "pairwise")
  expect_identical(lh$method, "haebara_pw")
})

test_that("ref re-expresses the trend in the reference metric", {
  ip <- toy_ipars()
  l1 <- link_chain(ip, method = "mgm")
  l2 <- link_chain(ip, method = "mgm", ref = 2)
  r <- match(2, l1$trend$group)
  expect_equal(l2$trend$mu[r], 0)
  expect_equal(l2$trend$sigma[r], 1)
  expect_equal(l2$trend$mu,
               (l1$trend$mu - l1$trend$mu[r]) / l1$trend$sigma[r])
  expect_equal(l2$trend$sigma, l1$trend$sigma / l1$trend$sigma[r])
  # identified parameters without DIF give the exact trend in the group-2 metric
  expect_equal(l2$trend$mu[1], (0 - 0.3) / 1.1, tolerance = 1e-8)

  lj <- link_joint(ip, method = "phl", ref = 3)
  expect_equal(lj$trend$mu[3], 0)
  expect_equal(lj$trend$sigma[3], 1)

  expect_error(link_chain(ip, method = "mgm", ref = 9), "must be one of")
})

test_that("ref works with the refit jackknife", {
  ip <- toy_ipars()
  ip$b[ip$group == 2 & ip$item == "I01"] <-
    ip$b[ip$group == 2 & ip$item == "I01"] + 0.5
  lk <- link_chain(ip, method = "mgm", ref = 2)
  le <- linking_error(lk, method = "jackknife")$le
  expect_equal(le$le_mu[le$group == 2], 0)
  expect_true(le$le_mu[le$group == 1] > 0)
  # ref = first group matches the default for every method
  lk1 <- link_chain(ip, method = "mgm", ref = 1)
  lk0 <- link_chain(ip, method = "mgm")
  expect_equal(linking_error(lk1, method = "sandwich_osw")$le,
               linking_error(lk0, method = "sandwich_osw")$le)
})

test_that("derivative-based LE methods support a non-default ref", {
  ip <- toy_ipars()
  # two DIF items so that all residuals are nonzero
  ip$b[ip$group == 2 & ip$item == "I01"] <-
    ip$b[ip$group == 2 & ip$item == "I01"] + 0.5
  ip$b[ip$group == 3 & ip$item == "I05"] <-
    ip$b[ip$group == 3 & ip$item == "I05"] - 0.4

  for (bauart in list(list(fn = link_chain, args = list(method = "mgm")),
                      list(fn = link_joint, args = list(method = "phl")))) {
    lk <- do.call(bauart$fn, c(list(x = ip, ref = 2), bauart$args))
    # the numeric and the analytic rebasing path must agree
    osw <- linking_error(lk, method = "sandwich_osw")$le
    osw_a <- linking_error(lk, method = "sandwich_osw_analytic")$le
    expect_equal(osw$le_mu, osw_a$le_mu, tolerance = 1e-5)
    expect_equal(osw$le_sigma, osw_a$le_sigma, tolerance = 1e-5)
    # reference group error free, group 1 carries its own error
    expect_equal(osw$le_mu[osw$group == 2], 0)
    expect_true(osw$le_mu[osw$group == 1] > 0)
    # the AJK approximates the refit jackknife in the ref metric too
    ajk <- linking_error(lk, method = "ajk")$le
    jk <- linking_error(lk, method = "jackknife")$le
    expect_equal(ajk$le_mu, jk$le_mu, tolerance = 0.15)
    # esw runs
    expect_silent(linking_error(lk, method = "sandwich_esw"))
    # bc variants give the full table with te^2 = se^2 + le^2
    V <- diag(0.001, 2 * nrow(lk$ipars))
    bc <- linking_error(lk, method = "sandwich_osw_bc", vcov = V)$le
    expect_equal(bc$te_mu, sqrt(bc$se_mu^2 + bc$le_mu^2))
    expect_equal(bc$le_mu[bc$group == 2], 0)
    jbc <- linking_error(lk, method = "jackknife_bc", vcov = V)$le
    expect_equal(jbc$le_mu[jbc$group == 2], 0)
    expect_true(jbc$se_mu[jbc$group == 1] > 0)
  }

  # closed-form SL AJK rebases the delete-one target vectors
  lk_sl <- link_chain(ip, method = "sl", ref = 2)
  ajk_sl <- linking_error(lk_sl, method = "ajk")$le
  jk_sl <- linking_error(lk_sl, method = "jackknife")$le
  expect_equal(ajk_sl$le_mu[ajk_sl$group == 2], 0)
  expect_true(ajk_sl$le_mu[ajk_sl$group == 1] > 0)
  expect_equal(ajk_sl$le_mu, jk_sl$le_mu, tolerance = 0.2)
  # ref = 1 falls back to the default
  lk_sl1 <- link_chain(ip, method = "sl")
  expect_equal(linking_error(link_chain(ip, method = "sl", ref = 1),
                             method = "ajk")$le,
               linking_error(lk_sl1, method = "ajk")$le)
})

test_that("non-uniform DIF shifts discriminations, not difficulties", {
  set.seed(42)
  d <- sim_trend_data(n_groups = 3, N = 50, I = 10, overlap = 1,
                      dif = "random", dif_par = "a", dif_sd_a = 0.2)
  expect_true(all(attr(d, "dif_effects") == 0))
  expect_true(any(attr(d, "dif_effects_a") != 0))

  set.seed(42)
  db <- sim_trend_data(n_groups = 3, N = 50, I = 10, overlap = 1,
                       dif = "random", dif_par = "both")
  expect_true(any(attr(db, "dif_effects") != 0))
  expect_true(any(attr(db, "dif_effects_a") != 0))

  expect_error(
    sim_trend_data(n_groups = 2, N = 50, I = 10, dif = "random",
                   dif_par = "a", model = "1PL"),
    "2PL")
})

test_that("dependent person samples share ids and correlate", {
  set.seed(7)
  d <- sim_trend_data(n_groups = 3, N = 400, I = 12, overlap = 1,
                      mu = c(0, 0.2, 0.4), sigma = c(1, 1, 1),
                      person_overlap = 0.5, person_cor = 0.9)
  ids <- attr(d, "person_ids")
  expect_identical(ids[[1]], rownames(d[[1]]))
  expect_identical(ids[[2]][1:200], ids[[1]][1:200])
  expect_identical(ids[[3]][1:200], ids[[2]][1:200])
  expect_identical(length(unique(unlist(ids))), 400L + 200L + 200L)
  # shared persons, sum scores of adjacent groups correlate
  s1 <- rowSums(d[[1]][1:200, ])
  s2 <- rowSums(d[[2]][1:200, ])
  expect_true(stats::cor(s1, s2) > 0.3)
  # independent default path, fresh persons per group
  set.seed(7)
  d0 <- sim_trend_data(n_groups = 2, N = 50, I = 10, overlap = 1)
  expect_identical(length(unique(unlist(attr(d0, "person_ids")))), 100L)
})

test_that("estimate_trend passes ref through", {
  skip_on_cran()
  set.seed(11)
  d <- sim_trend_data(n_groups = 2, N = 300, I = 10, overlap = 1,
                      mu = c(0, 0.3), sigma = c(1, 1.1))
  tr <- estimate_trend(d, calibration = "separate", approach = "chain",
                       link = "mgm", ref = 2)
  expect_equal(tr$trend$mu[tr$trend$group == 2], 0)
  expect_equal(tr$trend$sigma[tr$trend$group == 2], 1)
})
