test_that("detect_dif rmsd runs end to end on an em calibration", {
  set.seed(104)
  d <- sim_trend_data(n_groups = 2, N = 800, I = 15, dif = "unbalanced",
                      dif_pct = 0.3, dif_effect = 0.8)
  cal <- calibrate(d, engine = "em")
  dif <- detect_dif(cal, method = "rmsd", cutoff = 0.05)
  expect_true(length(dif$flagged) >= 1)
})

test_that("linking error with em vcov runs through jackknife_bc", {
  set.seed(105)
  d <- sim_trend_data(n_groups = 3, N = 800, I = 15, dif = "none")
  cal <- calibrate(d, engine = "em", keep_vcov = TRUE)
  lk <- link_chain(cal, method = "mgm")
  le <- linking_error(lk, method = "jackknife_bc", vcov = cal$vcov_ipars)
  expect_true(all(is.finite(le$le$se_mu[-1])))
  expect_true(all(le$le$se_mu[-1] > 0))
})

test_that("estimate_trend works with the em engine", {
  set.seed(106)
  d <- sim_trend_data(n_groups = 2, N = 500, I = 12, dif = "none")
  tr <- estimate_trend(d, link = "mgm", engine = "em")
  expect_s3_class(tr, "irtlink_trend")
  expect_equal(tr$trend$group, 1:2)
  expect_true(is.finite(tr$trend$mu[2]))
})
