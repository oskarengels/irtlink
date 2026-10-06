# tests/testthat/test-link-response-variants.R
# Symmetric and grid-metric variants of the Haebara / Stocking-Lord
# chain steps (symmetric and from-metric variants):
# `type = "symm"` sums both alignment directions; `theta_metric` names
# the group whose metric carries the quadrature grid and weights ("to" =
# the later group, the previous default; "from" = the earlier group, via
# swapped groups and inverted transformation).

response_var_ipars <- function(seed = 11) {
  ip <- make_identified_ipars(c(0, 0.35), c(1, 1.15), I = 12)
  set.seed(seed)
  sel <- ip$group == 2
  ip$b[sel] <- ip$b[sel] + stats::rnorm(12, sd = 0.25)
  ip$a[sel] <- ip$a[sel] * exp(stats::rnorm(12, sd = 0.15))
  ip
}

test_that("both variants recover an exact drift-free trend", {
  mu <- c(0, 0.3, 0.6); sigma <- c(1, 1.1, 1.2)
  ip <- make_identified_ipars(mu, sigma, I = 10)
  for (m in c("haebara", "sl")) {
    for (variant in list(list(type = "symm"),
                         list(theta_metric = "from"))) {
      lk <- do.call(link_chain, c(list(x = ip, method = m), variant))
      expect_equal(lk$trend$mu, mu, tolerance = 1e-3,
                   label = paste(m, names(variant)))
      expect_equal(lk$trend$sigma, sigma, tolerance = 1e-3,
                   label = paste(m, names(variant)))
    }
  }
})

test_that("the refit jackknife works for the new variants", {
  ip <- response_var_ipars()
  lk <- link_chain(ip, method = "haebara", type = "symm")
  jk <- linking_error(lk, method = "jackknife")$le
  expect_gt(jk$le_mu[2], 0)
  lk2 <- link_chain(ip, method = "sl", theta_metric = "from")
  jk2 <- linking_error(lk2, method = "jackknife")$le
  expect_gt(jk2$le_mu[2], 0)
})

test_that("the AJK covers the symm and from-metric response function variants", {
  ip <- response_var_ipars()
  for (variant in list(list(method = "haebara", type = "symm"),
                       list(method = "haebara", theta_metric = "from"),
                       list(method = "sl", type = "symm"),
                       list(method = "sl", theta_metric = "from"))) {
    lk <- do.call(link_chain, c(list(x = ip), variant))
    jk <- linking_error(lk, method = "jackknife")$le
    ajk <- linking_error(lk, method = "ajk")$le
    lab <- paste(unlist(variant), collapse = "/")
    # Haebara one-steps track the refit closely; the SL closed form
    # carries the documented Eq.-32 approximation error.
    tol <- if (variant$method == "haebara") 0.04 else 0.10
    expect_equal(ajk$le_mu[-1], jk$le_mu[-1], tolerance = tol,
                 label = paste(lab, "mu"))
    expect_equal(ajk$le_sigma[-1], jk$le_sigma[-1], tolerance = tol,
                 label = paste(lab, "sigma"))
  }
})

test_that("closed-form methods reject the response-function-only arguments", {
  ip <- make_identified_ipars(c(0, 0.3), c(1, 1.1), I = 6)
  expect_error(link_chain(ip, method = "mgm", type = "symm"), "haebara")
  expect_error(link_chain(ip, method = "mm", theta_metric = "from"),
               "haebara")
})
