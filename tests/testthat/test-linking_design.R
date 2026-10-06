test_that("linking_design builds item x group incidence matrix", {
  set.seed(1)
  d <- sim_trend_data(n_groups = 3, N = 10, I = 10, overlap = 0.4,
                      design = "all")
  ld <- linking_design(d)
  expect_s3_class(ld, "irtlink_design")
  expect_equal(ncol(ld$matrix), 3)
  expect_equal(nrow(ld$matrix), length(unique(unlist(lapply(d, colnames)))))
  expect_equal(ld$successive_common, c(4L, 4L))  # round(10 * 0.4)
})

test_that("linking_design detects zero overlap", {
  d <- list(
    data.frame(I1 = c(0, 1), I2 = c(1, 0)),
    data.frame(I3 = c(0, 1), I4 = c(1, 0))
  )
  ld <- linking_design(d)
  expect_equal(ld$successive_common, 0L)
})

test_that("linking_design validates input", {
  expect_error(linking_design(list(data.frame(I1 = 1))), "at least two")
  m <- matrix(0:1, 2, 2)  # no colnames
  expect_error(linking_design(list(m, m)), "column names")
})

test_that("print method runs without error", {
  set.seed(1)
  d <- sim_trend_data(n_groups = 2, N = 10, I = 10)
  expect_output(print(linking_design(d)), "2 groups")
})

test_that("linking_design handles a single-item pool", {
  d <- list(data.frame(I1 = c(0, 1)), data.frame(I1 = c(1, 0)))
  ld <- linking_design(d)
  expect_true(is.matrix(ld$matrix))
  expect_equal(dim(ld$matrix), c(1L, 2L))
  expect_equal(ld$successive_common, 1L)
})

test_that("print notes suppression for large pools", {
  d <- list(
    as.data.frame(matrix(0, 2, 61, dimnames = list(NULL, paste0("I", 1:61)))),
    as.data.frame(matrix(0, 2, 61, dimnames = list(NULL, paste0("I", 1:61))))
  )
  expect_output(print(linking_design(d)), "matrix suppressed")
})
