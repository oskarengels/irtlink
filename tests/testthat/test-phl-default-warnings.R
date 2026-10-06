# tests/testthat/test-phl-default-warnings.R

test_that("PHL defaults to inverse_admin item weights", {
  ip <- make_identified_ipars(c(0, 0.3, 0.6), c(1, 1.1, 1.2), I = 10)
  lk <- link_joint(ip, method = "phl")
  expect_identical(lk$control$item_weights, "inverse_admin")
  lk_explicit <- link_joint(ip, method = "phl",
                            item_weights = "inverse_admin")
  expect_equal(lk$trend, lk_explicit$trend)
  expect_identical(eval(formals(joint_phl)$item_weights)[1],
                   "inverse_admin")
  expect_identical(eval(formals(estimate_trend)$item_weights)[1],
                   "inverse_admin")
  expect_identical(le_phl_variant(list(control = list())),
                   "inverse_admin")
})

test_that("the jackknife warns when delete-one estimates are not finite", {
  ip <- make_identified_ipars(c(0, 0.3, 0.6), c(1, 1.1, 1.2), I = 6)
  lk <- link_chain(ip, method = "mgm")
  real_refit <- refit_link_without_item
  local_mocked_bindings(
    refit_link_without_item = function(link, item, ...) {
      out <- real_refit(link, item, ...)
      if (identical(item, "I01")) out$trend$mu[2] <- NaN
      out
    }
  )
  expect_warning(res <- linking_error(lk, method = "jackknife"),
                 "not finite")
  le <- res$le
  expect_true(all(is.finite(le$le_mu[le$group > 1])))
})

test_that("estimate_trend reports the independence fallback for *_bc", {
  set.seed(11)
  d <- sim_trend_data(n_groups = 2, N = 150, I = 8, overlap = 1,
                      dif = "none", mu = c(0, 0.3), sigma = c(1, 1))
  expect_message(
    estimate_trend(d, calibration = "separate", approach = "chain",
                   link = "mgm", le = "jackknife_bc"),
    "independent")
  expect_no_message(
    estimate_trend(d, calibration = "separate", approach = "chain",
                   link = "mgm", le = "jackknife"))
})
