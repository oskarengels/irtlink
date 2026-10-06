# tests/testthat/test-link-methods-plot.R

test_that("plot.irtlink_link draws drift panels without error", {
  ip <- make_identified_ipars(c(0, 0.3, 0.6), c(1, 1.1, 1.2))
  res <- link_chain(ip, method = "mgm")
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  expect_no_error(plot(res))
  expect_no_error(plot(res, pars = "a"))
  resj <- link_joint(ip, method = "phl")
  expect_no_error(plot(resj))
})

test_that("plot errors informatively on objects without ipars", {
  ip <- make_identified_ipars(c(0, 0.3), c(1, 1))
  res <- link_chain(ip, method = "mgm")
  res$ipars <- NULL
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  expect_error(plot(res), "no item parameters")
})

test_that("plot forwards graphical args including main/xlab/ylab", {
  ip <- make_identified_ipars(c(0, 0.3), c(1, 1.1))
  res <- link_chain(ip, method = "mgm")
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  expect_no_error(plot(res, col = "blue"))
  expect_no_error(plot(res, main = "custom"))
  expect_no_error(plot(res, xlab = "x", ylab = "y"))
})
