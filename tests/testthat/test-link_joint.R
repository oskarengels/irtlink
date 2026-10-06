# tests/testthat/test-link_joint.R

joint_recovery_check <- function(tol) {
  mu <- c(0, 0.3, 0.6); sigma <- c(1, 1.1, 1.2)
  ip <- make_identified_ipars(mu, sigma)
  for (m in names(tol)) {
    res <- suppressWarnings(link_joint(ip, method = m))
    expect_s3_class(res, "irtlink_link")
    expect_identical(res$approach, "joint")
    expect_null(res$steps)
    expect_equal(res$trend$mu, mu, tolerance = tol[[m]], label = paste("mu", m))
    expect_equal(res$trend$sigma, sigma, tolerance = tol[[m]],
                 label = paste("sigma", m))
  }
}

test_that("link_joint recovers a 3-group trend for all methods", {
  joint_recovery_check(c(phl = 1e-8, haebara = 1e-3, haberman = 1e-6,
                         sl = 1e-4))
})

test_that("link_joint orders groups correctly beyond 9 groups", {
  mu <- 0.1 * (0:10); sigma <- rep(1, 11)
  ip <- make_identified_ipars(mu, sigma, I = 10)
  res <- link_joint(ip, method = "phl")
  expect_equal(res$trend$mu, mu, tolerance = 1e-8)
})

test_that("link_joint rejects chain-only methods", {
  ip <- make_identified_ipars(c(0, 0.3), c(1, 1))
  expect_error(link_joint(ip, method = "mgm"), "should be one of")
})

test_that("anchor restricts items and methods print without steps", {
  ip <- make_identified_ipars(c(0, 0.3), c(1, 1.1), I = 10)
  res <- suppressWarnings(
    link_joint(ip, method = "haberman", anchor = sprintf("I%02d", 1:5))
  )
  expect_equal(sort(unique(res$ipars$item)), sprintf("I%02d", 1:5))
  expect_output(print(res), "joint linking \\(haberman\\)")
  expect_output(summary(res), "joint")
  expect_identical(coef(res), res$trend)
})

test_that("link_chain objects now carry their ipars", {
  ip <- make_identified_ipars(c(0, 0.3), c(1, 1.1))
  res <- link_chain(ip, method = "mgm")
  expect_identical(res$ipars, as_ipars(ip))
})

test_that("link_joint validates pow like link_chain", {
  ip <- make_identified_ipars(c(0, 0.3), c(1, 1.1))
  expect_error(link_joint(ip, method = "haberman", pow = 3),
               "`pow` must be one of")
  expect_error(link_joint(ip, method = "haebara", pow = 3),
               "`pow` must be one of")
})

test_that("link_joint warns when anchor filtering drops a whole group", {
  # group 2 keeps only items I01-I05 after the drop below; anchoring on I06-I10
  # then removes group 2 entirely
  ip2 <- make_identified_ipars(c(0, 0.3, 0.6), c(1, 1.1, 1.2), I = 10)
  ip2 <- ip2[!(ip2$group == 2 & ip2$item %in% sprintf("I%02d", 6:10)), ]
  expect_warning(
    link_joint(ip2, method = "phl", anchor = sprintf("I%02d", 6:10)),
    "removed all items from group"
  )
})

test_that("restricted joint recovers the trend and matches subset runs", {
  mu <- c(0, 0.3, 0.6, 0.9); sigma <- c(1, 1.1, 1.2, 1.3)
  ip <- make_identified_ipars(mu, sigma)
  res <- link_joint(ip, method = "phl", restricted = TRUE)
  expect_identical(res$approach, "joint_restricted")
  expect_length(res$converged, 4)
  expect_equal(res$trend$mu, mu, tolerance = 1e-8)
  expect_equal(res$trend$sigma, sigma, tolerance = 1e-8)
  # group-2 estimate equals a plain 2-group joint over groups 1..2
  res2 <- link_joint(ip[ip$group <= 2, ], method = "phl")
  expect_equal(res$trend$mu[2], res2$trend$mu[2], tolerance = 1e-10)
})

test_that("link_joint rejects non-positive discriminations for all backends", {
  ip <- as_ipars(make_identified_ipars(c(0, 0.3, 0.6), c(1, 1.1, 1.2), I = 10))
  ip$a[1] <- -0.5
  for (m in c("haberman", "haebara", "sl", "phl")) {
    expect_error(suppressWarnings(link_joint(ip, method = m)),
                 "Non-positive discrimination",
                 label = m)
  }
})

test_that("link_joint exposes eps override for haberman", {
  ip <- as_ipars(make_identified_ipars(c(0, 0.3), c(1, 1.1)))
  res <- suppressWarnings(link_joint(ip, method = "haberman", eps = 0.01))
  expect_equal(res$trend$mu[2], 0.3, tolerance = 1e-3)
})
