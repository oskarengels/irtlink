# tests/testthat/test-calibrate.R
# TAM calibrations are moderately slow; keep N and I small.

make_test_data <- function(n_groups = 2, N = 300, I = 10, overlap = 0.5) {
  set.seed(42)
  sim_trend_data(n_groups = n_groups, N = N, I = I, overlap = overlap,
                 design = "all")
}

test_that("calibrate returns a valid irtlink_calib object", {
  d <- make_test_data()
  cal <- calibrate(d, model = "2PL")
  expect_s3_class(cal, "irtlink_calib")
  expect_identical(cal$model, "2PL")
  expect_identical(cal$calibration, "separate")
  expect_identical(cal$engine, "em")          # the internal engine default
  expect_s3_class(cal$models[[1]], "irtlink_em")
  expect_equal(cal$n_groups, 2)
  expect_equal(cal$N, c(300L, 300L))
  # one ipars row per administered item per group
  expect_equal(nrow(cal$ipars), 20)
  expect_identical(names(cal$ipars), c("group", "item", "a", "b", "c"))
  expect_length(cal$models, 2)
  # as_ipars() unwraps the calib object
  expect_identical(as_ipars(cal), cal$ipars)
  # ipars must satisfy the canonical contract (validated via as_ipars)
  expect_identical(as_ipars(cal$ipars), cal$ipars)
})

test_that("calibrate recovers item difficulties at group 1", {
  d <- make_test_data()
  truth <- attr(d, "true_ipars")
  cal <- calibrate(d, model = "2PL")
  w1 <- cal$ipars[cal$ipars$group == 1, ]
  truth_w1 <- truth[match(w1$item, truth$item), ]
  expect_gt(cor(w1$b, truth_w1$b), 0.9)
})

test_that("keep_models = FALSE drops the fitted model objects", {
  d <- make_test_data()
  cal <- calibrate(d, keep_models = FALSE)
  expect_null(cal$models)
})

test_that("separate calibration can store item-parameter vcov matrices", {
  skip_on_cran()
  set.seed(421)
  d <- sim_trend_data(n_groups = 2, N = 250, I = 6, overlap = 1,
                      model = "2PL")
  cal <- calibrate(d, model = "2PL", keep_vcov = TRUE)

  expect_length(cal$vcov_ipars, 2)
  for (t in seq_along(cal$vcov_ipars)) {
    V <- cal$vcov_ipars[[t]]
    n_items <- ncol(d[[t]])
    expect_equal(dim(V), c(2 * n_items, 2 * n_items))
    expect_true(all(paste0(colnames(d[[t]]), "_a") %in% rownames(V)))
    expect_true(all(paste0(colnames(d[[t]]), "_b") %in% rownames(V)))
    expect_lt(max(abs(V - t(V))), 1e-8)
  }
})

test_that("calibrate validates inputs", {
  d <- make_test_data()
  bad <- d
  bad[[1]][1, 1] <- 2
  expect_error(calibrate(bad), "non-dichotomous")
  no_overlap <- list(
    data.frame(I1 = rbinom(50, 1, 0.5), I2 = rbinom(50, 1, 0.5)),
    data.frame(I3 = rbinom(50, 1, 0.5), I4 = rbinom(50, 1, 0.5))
  )
  expect_error(calibrate(no_overlap), "no common items")
})

test_that("calibrate warns on thin overlap (< 5 common items)", {
  set.seed(7)
  d <- sim_trend_data(n_groups = 2, N = 200, I = 10, overlap = 0.3)
  expect_warning(calibrate(d), "Fewer than 5 common items")
})

test_that("1PL calibration fixes discriminations at 1", {
  set.seed(42)
  d <- sim_trend_data(n_groups = 2, N = 300, I = 10, overlap = 0.5,
                      model = "1PL")
  cal <- calibrate(d, model = "1PL")
  expect_true(all(cal$ipars$a == 1))
  truth <- attr(d, "true_ipars")
  w1 <- cal$ipars[cal$ipars$group == 1, ]
  expect_gt(cor(w1$b, truth$b[match(w1$item, truth$item)]), 0.9)
})
