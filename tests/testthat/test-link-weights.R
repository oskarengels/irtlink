# tests/testthat/test-link-weights.R

test_that("make_link_weights builds the quadrature weight schemes", {
  w <- make_link_weights("uniform")
  expect_length(w$theta, 101)
  expect_equal(range(w$theta), c(-6, 6))
  expect_true(all(w$wgt == 1))

  w05 <- make_link_weights("normal_0.5")
  expect_equal(w05$wgt, dnorm(w05$theta, mean = 0, sd = 0.5))
  expect_equal(which.max(w05$wgt), 51L)  # peak at theta = 0
})

test_that("make_link_weights validates the weights label", {
  expect_error(make_link_weights("triangular"), "should be one of")
})
