# tests/testthat/test-detect_dif-iterative.R

test_that("iterative purification removes drifted items and records history", {
  cal <- make_drift_calib(drift_items = c("I05", "I12"), drift = 1.2)
  res <- detect_dif(cal, rmsd_method = "parameter", cutoff = 0.05,
                    procedure = "iterative_forward")
  expect_true(all(c("I05", "I12") %in% res$flagged))
  expect_false(any(c("I05", "I12") %in% res$anchor))
  expect_s3_class(res, "irtlink_dif")
  expect_true(is.list(res$history))
  expect_gte(length(res$history), 1)
})

test_that("iterative stops at min_anchor and never over-flags", {
  cal <- make_drift_calib(drift_items = sprintf("I%02d", 1:8), drift = 1.5,
                          I = 10)
  res <- detect_dif(cal, rmsd_method = "parameter", cutoff = 0.05,
                    procedure = "iterative_forward", min_anchor = 3)
  # only items present in >= 2 groups are linkable; anchors never below 3
  expect_gte(length(res$anchor), 3)
})

test_that("iterative with no drift flags nothing in one round", {
  cal <- make_drift_calib(drift_items = character(0))
  res <- detect_dif(cal, rmsd_method = "parameter", cutoff = 0.05,
                    procedure = "iterative_forward")
  expect_length(res$flagged, 0)
  expect_length(res$history, 1)
})

test_that("iterative warns and truncates a cutoff vector", {
  cal <- make_drift_calib(drift_items = "I05", drift = 1.0)
  expect_warning(
    detect_dif(cal, rmsd_method = "parameter", cutoff = c(0.03, 0.08),
               procedure = "iterative_forward"),
    "one_step"
  )
})

test_that("iterative warns when max_iter is reached before convergence", {
  cal <- make_drift_calib(drift_items = sprintf("I%02d", 1:6), drift = 1.3,
                          I = 20)
  expect_warning(
    detect_dif(cal, rmsd_method = "parameter", cutoff = 0.05,
               procedure = "iterative_forward", max_iter = 1),
    "max_iter"
  )
})

test_that("iterative result reports convergence", {
  cal_clean <- make_drift_calib(drift_items = "I05", drift = 1.0)
  res <- detect_dif(cal_clean, rmsd_method = "parameter", cutoff = 0.05,
                    procedure = "iterative_forward")
  expect_true(res$converged)
})

test_that("both + iterative_forward warns that the other estimator's table is NULL", {
  cal <- make_drift_calib(drift_items = "I05", drift = 1.0)
  expect_warning(
    detect_dif(cal, rmsd_method = "both", flag_on = "parameter",
               procedure = "iterative_forward"),
    "the other estimator's table is NULL"
  )
})
