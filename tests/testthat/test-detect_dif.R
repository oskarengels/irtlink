# tests/testthat/test-detect_dif.R

test_that("detect_dif (parameter, one_step) flags the drifted item", {
  cal <- make_drift_calib(drift_items = "I05", drift = 1.0)
  res <- detect_dif(cal, method = "rmsd", rmsd_method = "parameter",
                    cutoff = 0.05, procedure = "one_step")
  expect_s3_class(res, "irtlink_dif")
  expect_true("I05" %in% res$flagged)
  expect_false("I05" %in% res$anchor)
  expect_identical(res$rmsd_method, "parameter")
})

test_that("one_step accepts a cutoff vector and stores per-cutoff flags", {
  cal <- make_drift_calib(drift_items = "I05", drift = 1.0)
  res <- detect_dif(cal, rmsd_method = "parameter",
                    cutoff = c(0.03, 0.08), procedure = "one_step")
  expect_named(res$flagged_by_cutoff, c("0.03", "0.08"))
  # the primary flagged set uses the first cutoff
  expect_identical(res$flagged, res$flagged_by_cutoff[["0.03"]])
})

test_that("data-driven cutoff path runs on parameter RMSD", {
  cal <- make_drift_calib(drift_items = "I05", drift = 1.2)
  res <- detect_dif(cal, rmsd_method = "parameter", cutoff_type = "data_driven",
                    tau = 2.7, procedure = "one_step")
  expect_true("I05" %in% res$flagged)
})

test_that("lrt and data+iterative both require keep_models = TRUE", {
  cal <- make_drift_calib(drift_items = character(0), mu = c(0, 0.3),
                          sigma = c(1, 1.1))
  # LRT is now implemented; with models = NULL it errors asking for keep_models
  expect_error(detect_dif(cal, method = "lrt"), "fitted models")
  # data + iterative: same requirement
  expect_error(
    detect_dif(cal, rmsd_method = "data", procedure = "iterative_forward"),
    "fitted models"
  )
})

test_that("flag_on is ignored with a warning when rmsd_method is not both", {
  cal <- make_drift_calib(drift_items = "I05", drift = 1.0)
  expect_warning(
    detect_dif(cal, rmsd_method = "parameter", flag_on = "data"),
    "ignored"
  )
})

test_that("detect_dif rejects a non-calib input", {
  expect_error(detect_dif(data.frame(a = 1)), "irtlink_calib")
})
