# tests/testthat/test-calib-methods.R

make_calib_fixture <- function() {
  # Hand-built object: avoids slow TAM fits for method tests.
  ipars <- data.frame(
    group = rep(1:2, each = 3),
    item = rep(c("I1", "I2", "I3"), 2),
    a = c(1.0, 1.2, 0.8, 1.1, 1.3, 0.9),
    b = c(-0.5, 0.0, 0.5, -0.4, 0.1, 0.6)
  )
  structure(
    list(ipars = ipars, model = "2PL", calibration = "separate",
         engine = "em", n_groups = 2L, N = c(100L, 100L),
         converged = c(TRUE, TRUE), models = NULL, design = NULL,
         call = NULL),
    class = "irtlink_calib"
  )
}

test_that("print shows model, groups, and convergence", {
  cal <- make_calib_fixture()
  expect_output(print(cal), "2PL")
  expect_output(print(cal), "Groups: 2")
  cal$converged <- c(TRUE, FALSE)
  expect_output(print(cal), "did not converge.*2")
  cal$converged <- c(TRUE, NA)
  expect_output(print(cal), "did not converge.*2")
})

test_that("coef returns the ipars data.frame", {
  cal <- make_calib_fixture()
  expect_identical(coef(cal), cal$ipars)
})

test_that("summary prints parameter ranges per group", {
  cal <- make_calib_fixture()
  expect_output(summary(cal), "group")
  expect_output(summary(cal), "a")
  expect_output(summary(cal), "0\\.8")
})
