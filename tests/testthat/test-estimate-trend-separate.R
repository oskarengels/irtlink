test_that("estimate_trend separate+chain recovers the trend", {
  skip_on_cran()
  set.seed(4101)
  d <- sim_trend_data(n_groups = 3, N = 2000, I = 15, dif = "none",
                      mu = c(0, 0.4, 0.8), sigma = c(1, 1.1, 1.2))
  tr <- estimate_trend(d, calibration = "separate", approach = "chain",
                       link = "mgm")
  expect_s3_class(tr, "irtlink_trend")
  expect_named(tr$trend, c("group", "mu", "sigma", "se", "le", "te"))
  expect_equal(tr$trend$mu, c(0, 0.4, 0.8), tolerance = 0.12)
  expect_equal(tr$trend$sigma, c(1, 1.1, 1.2), tolerance = 0.12)
  expect_true(all(is.na(tr$trend$le)))
  expect_s3_class(tr$calib, "irtlink_calib")
  expect_s3_class(tr$link_obj, "irtlink_link")
})

test_that("estimate_trend separate+joint and joint_restricted run", {
  skip_on_cran()
  set.seed(4102)
  d <- sim_trend_data(n_groups = 3, N = 2000, I = 15, dif = "none",
                      mu = c(0, 0.3, 0.6), sigma = c(1, 1, 1))
  trj <- estimate_trend(d, "separate", "joint", link = "haberman")
  trr <- estimate_trend(d, "separate", "joint_restricted", link = "haberman")
  expect_equal(trj$trend$mu, c(0, 0.3, 0.6), tolerance = 0.12)
  expect_equal(trr$trend$mu, c(0, 0.3, 0.6), tolerance = 0.12)
  expect_identical(trj$approach, "joint")
  expect_identical(trr$approach, "joint_restricted")
})

test_that("mm/mgm are rejected for joint approaches", {
  d <- sim_trend_data(n_groups = 2, N = 300, I = 10, dif = "none")
  expect_error(estimate_trend(d, "separate", "joint", link = "mgm"), "chain")
})

test_that("estimate_trend separate+chain can attach jackknife linking error", {
  skip_on_cran()
  set.seed(4103)
  d <- sim_trend_data(n_groups = 2, N = 800, I = 12, dif = "none",
                      mu = c(0, 0.3), sigma = c(1, 1.1))
  tr <- estimate_trend(d, "separate", "chain", link = "mgm",
                       le = "jackknife")
  expect_s3_class(tr, "irtlink_trend")
  expect_false(is.null(tr$le_obj))
  expect_true(any(!is.na(tr$trend$le)))
  expect_identical(tr$le_obj$method, rep("jackknife", nrow(tr$le_obj)))
})

test_that("estimate_trend separate+chain can attach AJK linking error", {
  skip_on_cran()
  set.seed(4105)
  d <- sim_trend_data(n_groups = 2, N = 800, I = 12, dif = "none",
                      mu = c(0, 0.3), sigma = c(1, 1.1))
  tr <- estimate_trend(d, "separate", "chain", link = "mgm",
                       le = "ajk")
  expect_s3_class(tr, "irtlink_trend")
  expect_false(is.null(tr$le_obj))
  expect_true(any(!is.na(tr$trend$le)))
  expect_identical(tr$le_obj$method, rep("ajk", nrow(tr$le_obj)))
})

test_that("estimate_trend separate+chain can attach ESW linking error", {
  skip_on_cran()
  set.seed(4104)
  d <- sim_trend_data(n_groups = 2, N = 800, I = 12, dif = "none",
                      mu = c(0, 0.3), sigma = c(1, 1.1))
  tr <- estimate_trend(d, "separate", "chain", link = "mgm",
                       le = "sandwich_esw")
  expect_s3_class(tr, "irtlink_trend")
  expect_false(is.null(tr$le_obj))
  expect_true(any(!is.na(tr$trend$le)))
  expect_identical(tr$le_obj$method, rep("sandwich_esw", nrow(tr$le_obj)))
})

test_that("estimate_trend separate+chain can attach analytic OSW linking error", {
  skip_on_cran()
  set.seed(4106)
  d <- sim_trend_data(n_groups = 2, N = 800, I = 12, dif = "none",
                      mu = c(0, 0.3), sigma = c(1, 1.1))
  tr <- estimate_trend(d, "separate", "chain", link = "mgm",
                       le = "sandwich_osw_analytic")
  expect_s3_class(tr, "irtlink_trend")
  expect_false(is.null(tr$le_obj))
  expect_true(any(!is.na(tr$trend$le)))
  expect_identical(tr$le_obj$method,
                   rep("sandwich_osw_analytic", nrow(tr$le_obj)))
})

test_that("estimate_trend auto-uses the calibration vcov for bias-corrected AJK", {
  skip_on_cran()
  set.seed(4107)
  d <- sim_trend_data(n_groups = 2, N = 250, I = 6, overlap = 1,
                      dif = "none", mu = c(0, 0.3), sigma = c(1, 1.1))
  tr <- estimate_trend(d, "separate", "chain", link = "mgm",
                       le = "ajk_bc")

  expect_false(is.null(tr$calib$vcov_ipars))
  expect_true(all(c("se", "le", "te", "le_bc", "te_bc") %in%
                    names(tr$trend)))
  expect_true(any(!is.na(tr$trend$se)))
  expect_true(any(!is.na(tr$trend$te)))
  expect_identical(tr$le_obj$method, rep("ajk_bc", nrow(tr$le_obj)))
})

test_that("estimate_trend routes the pairwise joint methods and weights", {
  skip_on_cran()
  set.seed(4109)
  d <- sim_trend_data(n_groups = 2, N = 300, I = 10, overlap = 1,
                      dif = "none", mu = c(0, 0.3), sigma = c(1, 1.1))
  tr <- estimate_trend(d, "separate", "joint", link = "haebara",
                       variant = "pairwise",
                       weights = "uniform", le = "jackknife")
  expect_identical(tr$link_obj$method, "haebara_pw")
  expect_identical(tr$link_obj$control$weights, "uniform")
  expect_true(any(!is.na(tr$trend$le)))
  expect_error(estimate_trend(d, "separate", "chain", link = "sl",
                              variant = "pairwise"),
               "pairwise")
})

test_that("estimate_trend supports the bias-corrected sandwich linking error", {
  skip_on_cran()
  set.seed(4108)
  d <- sim_trend_data(n_groups = 2, N = 250, I = 6, overlap = 1,
                      dif = "none", mu = c(0, 0.3), sigma = c(1, 1.1))
  # em engine (default) stores vcov_ipars automatically for the bc methods
  tr <- estimate_trend(d, "separate", "chain", link = "mgm",
                       le = "sandwich_osw_bc")
  expect_identical(tr$le_obj$method,
                   rep("sandwich_osw_bc", nrow(tr$le_obj)))
  expect_true(all(c("se", "le", "te", "le_bc", "te_bc") %in%
                    names(tr$trend)))
  expect_true(any(!is.na(tr$trend$se)))
})
