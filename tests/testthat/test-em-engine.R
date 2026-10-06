test_that("calibrate dispatches to the em engine for separate", {
  set.seed(41)
  d <- sim_trend_data(n_groups = 2, N = 400, I = 15, dif = "none")
  cal <- calibrate(d, engine = "em")
  expect_s3_class(cal, "irtlink_calib")
  expect_equal(cal$engine, "em")
  expect_equal(sort(unique(cal$ipars$group)), 1:2)
  expect_named(cal$ipars, c("group", "item", "a", "b", "c"))
  expect_length(cal$models, 2L)
  expect_s3_class(cal$models[[1]], "irtlink_em")
  expect_true(all(cal$converged))
})

test_that("control options reach the em engine through ...", {
  set.seed(42)
  d <- sim_trend_data(n_groups = 2, N = 300, I = 15, dif = "none")
  cal <- calibrate(d, engine = "em", control = list(theta_nodes = 61L))
  expect_length(cal$models[[1]]$theta, 61L)
})

test_that("raw responses are recoverable from em fits", {
  set.seed(43)
  d <- sim_trend_data(n_groups = 2, N = 300, I = 15, dif = "none")
  cal <- calibrate(d, engine = "em")
  rpw <- responses_per_group(cal)
  expect_length(rpw, 2L)
  expect_equal(sort(colnames(rpw[[1]])),
               sort(colnames(as.data.frame(d[[1]]))))
  expect_equal(nrow(rpw[[1]]), nrow(as.data.frame(d[[1]])))
})

test_that("the em engine is the calibrate default", {
  set.seed(48)
  d <- sim_trend_data(n_groups = 2, N = 200, I = 15, dif = "none")
  cal <- calibrate(d)
  expect_equal(cal$engine, "em")
  expect_s3_class(cal$models[[1]], "irtlink_em")
})

test_that("calibrate supports the 3PL on the em engine only", {
  set.seed(49)
  d <- sim_trend_data(n_groups = 2, N = 400, I = 15, dif = "none")
  cal <- calibrate(d, model = "3PL")
  expect_equal(cal$model, "3PL")
  expect_named(cal$ipars, c("group", "item", "a", "b", "c"))
  expect_length(cal$models[[1]]$c, 15L)
  expect_error(calibrate(d, model = "3PL", calibration = "concurrent"),
               "3PL")
  expect_error(person_scores(cal), "3PL")
})

test_that("calibrate concurrent dispatches to the em engine", {
  set.seed(45)
  d <- sim_trend_data(n_groups = 3, N = 400, I = 15, dif = "none")
  cal <- calibrate(d, calibration = "concurrent", engine = "em")
  expect_equal(cal$engine, "em")
  expect_equal(cal$trend$group, 1:3)
  expect_equal(cal$trend$mu[1], 0)
  expect_equal(cal$trend$sigma[1], 1)
  expect_named(cal$ipars, c("group", "item", "a", "b", "c"))
  expect_length(cal$models, 1L)
  expect_s3_class(cal$models[[1]], "irtlink_em")
  # invariant parameters replicated across the groups each item appears in
  ip <- cal$ipars
  for (it in unique(ip$item)) {
    expect_length(unique(ip$a[ip$item == it]), 1L)
    expect_length(unique(ip$b[ip$item == it]), 1L)
  }
})

test_that("calibrate concurrent em supports free_items", {
  set.seed(46)
  d <- sim_trend_data(n_groups = 2, N = 400, I = 15, dif = "none")
  it <- colnames(as.data.frame(d[[1]]))[1]
  cal <- calibrate(d, calibration = "concurrent", engine = "em",
                   free_items = it)
  expect_equal(cal$free_items, it)
  expect_true(paste0(it, "__G2_g") %in%
                names(cal$models[[1]]$g_estimates))
})

test_that("extract_em_ipars warns on bound-stuck discrimination", {
  fit <- list(item = c("i1", "i2"), a = c(0.1, 1), b = c(0, 0))
  expect_warning(extract_em_ipars(fit, group = 1L), "discrimination bound")
})
