# tests/testthat/test-joint-phl.R

test_that("joint_phl recovers a 3-group trend exactly (closed form)", {
  mu <- c(0, 0.3, 0.6); sigma <- c(1, 1.1, 1.2)
  ip <- as_ipars(make_identified_ipars(mu, sigma))
  for (v in c("uniform", "inverse_admin")) {
    est <- joint_phl(ip, item_weights = v)
    expect_equal(est$mu, mu, tolerance = 1e-8)
    expect_equal(est$sigma, sigma, tolerance = 1e-8)
  }
})

test_that("the two item weightings differ under unbalanced item presence", {
  set.seed(8)
  ip <- as_ipars(make_identified_ipars(c(0, 0.3, 0.6), c(1, 1, 1), I = 12))
  # remove some items from single groups -> G_i varies -> omega differs
  drop <- (ip$group == 1 & ip$item %in% c("I01", "I02")) |
          (ip$group == 3 & ip$item %in% c("I03", "I04"))
  ip <- ip[!drop, ]
  ip$b <- ip$b + rnorm(nrow(ip), sd = 0.15)
  e1 <- joint_phl(ip, "uniform"); e2 <- joint_phl(ip, "inverse_admin")
  expect_false(isTRUE(all.equal(e1$mu, e2$mu, tolerance = 1e-10)))
})

test_that("joint_phl rejects non-positive discriminations", {
  ip <- as_ipars(make_identified_ipars(c(0, 0.3), c(1, 1.1), I = 10))
  ip$a[3] <- -0.5
  expect_error(joint_phl(ip, "uniform"), "Non-positive discrimination")
})

test_that("joint_phl reports disconnected group graphs informatively", {
  # groups 1-2 and 3-4 are linked internally but no item connects the groups
  ip <- data.frame(
    group = c(1, 1, 2, 2, 3, 3, 4, 4),
    item = c("A", "B", "A", "B", "C", "D", "C", "D"),
    a = rep(1, 8),
    b = c(0, 0.5, -0.3, 0.2, 0, 0.5, -0.3, 0.2)
  )
  expect_error(joint_phl(as_ipars(ip), "uniform"), "disconnected")
})
