# tests/testthat/test-link_chain.R

test_that("link_chain recovers a 4-group trend exactly (identified pars)", {
  mu <- c(0, 0.3, 0.6, 0.9); sigma <- c(1, 1.1, 1.2, 1.3)
  ip <- make_identified_ipars(mu, sigma)
  for (m in c("mm", "mgm")) {
    res <- link_chain(ip, method = m)
    expect_s3_class(res, "irtlink_link")
    expect_equal(res$trend$mu, mu, tolerance = 1e-8)
    expect_equal(res$trend$sigma, sigma, tolerance = 1e-8)
  }
})

test_that("link_chain haberman recovers a 4-group trend", {
  mu <- c(0, 0.3, 0.6, 0.9); sigma <- c(1, 1.1, 1.2, 1.3)
  ip <- make_identified_ipars(mu, sigma)
  res_hab <- suppressWarnings(link_chain(ip, method = "haberman"))
  expect_equal(res_hab$trend$mu, mu, tolerance = 1e-4)
  expect_equal(res_hab$trend$sigma, sigma, tolerance = 1e-4)
})

test_that("group 1 is fixed at mu = 0, sigma = 1", {
  ip <- make_identified_ipars(c(0, 0.5), c(1, 1.2))
  res <- link_chain(ip, method = "mgm")
  expect_equal(res$trend$mu[1], 0)
  expect_equal(res$trend$sigma[1], 1)
})

test_that("steps table records each successive pair", {
  ip <- make_identified_ipars(c(0, 0.3, 0.6), c(1, 1, 1))
  res <- link_chain(ip, method = "mgm")
  expect_equal(res$steps$from, c(1, 2))
  expect_equal(res$steps$to, c(2, 3))
  expect_equal(res$steps$n_common, c(20L, 20L))
  expect_true(all(res$steps$converged))
})

test_that("anchor argument restricts the linking items", {
  ip <- make_identified_ipars(c(0, 0.3), c(1, 1.1))
  res <- link_chain(ip, method = "mgm", anchor = c("I01", "I02", "I03"))
  expect_equal(res$steps$n_common, 3L)
  expect_equal(res$trend$mu[2], 0.3, tolerance = 1e-8)
  expect_error(link_chain(ip, method = "mgm", anchor = "NOPE"),
               "No anchor items")
})

test_that("link_chain validates input", {
  ip <- make_identified_ipars(c(0), c(1))
  expect_error(link_chain(ip), "consecutive integers|at least two")
  ip2 <- make_identified_ipars(c(0, 0.3), c(1, 1))
  expect_error(link_chain(ip2, method = "phl"), "should be one of")
})

test_that("print, summary, and coef methods work", {
  ip <- make_identified_ipars(c(0, 0.3), c(1, 1.1))
  res <- link_chain(ip, method = "mgm")
  expect_output(print(res), "chain linking \\(mgm\\)")
  expect_output(summary(res), "Step details")
  expect_identical(coef(res), res$trend)
})

test_that("partially missing anchor items produce a message", {
  ip <- make_identified_ipars(c(0, 0.3), c(1, 1.1), I = 10)
  expect_message(
    link_chain(ip, method = "mgm", anchor = c("I01", "I02", "NOPE")),
    "2 of 3 anchor item"
  )
})

test_that("print warns about non-converged steps", {
  ip <- make_identified_ipars(c(0, 0.3), c(1, 1.1))
  res <- link_chain(ip, method = "mgm")
  res$steps$converged <- FALSE
  expect_output(print(res), "WARNING: step\\(s\\) did not converge")
})

test_that("link_chain supports haebara and sl with weights", {
  mu <- c(0, 0.3, 0.6); sigma <- c(1, 1.1, 1.2)
  ip <- make_identified_ipars(mu, sigma)
  for (m in c("haebara", "sl")) {
    res <- link_chain(ip, method = m, weights = "normal_1")
    expect_equal(res$trend$mu, mu, tolerance = 1e-3)
    expect_equal(res$trend$sigma, sigma, tolerance = 1e-3)
  }
})
