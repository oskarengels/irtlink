# guessing parameters in the response function linking criteria

ip_3pl <- function(cval, I = 12) {
  b_ref <- seq(-1.5, 1.5, length.out = I)
  a_ref <- rep(c(0.8, 1.2, 1.0), length.out = I)
  data.frame(
    group = rep(1:2, each = I),
    item = rep(sprintf("I%02d", seq_len(I)), 2),
    a = c(a_ref, a_ref * 1.1),
    b = c(b_ref, (b_ref - 0.3) / 1.1),
    c = cval
  )
}

test_that("haebara and sl recover an exact 3PL transformation", {
  ip <- ip_3pl(0.2)
  for (m in c("haebara", "sl")) {
    lk <- link_groups(ip, approach = "chain", method = m)
    expect_equal(lk$trend$mu[2], 0.3, tolerance = 1e-6)
    expect_equal(lk$trend$sigma[2], 1.1, tolerance = 1e-6)
  }
  lk <- link_groups(ip, approach = "chain", method = "haebara",
                    type = "symm")
  expect_equal(lk$trend$mu[2], 0.3, tolerance = 1e-6)
  expect_equal(lk$trend$sigma[2], 1.1, tolerance = 1e-6)
})

test_that("c = 0 reproduces the 2PL response function link exactly", {
  ip <- ip_3pl(0)
  ip$b[ip$group == 2][1:3] <- ip$b[ip$group == 2][1:3] + 0.4
  ip_noc <- ip[c("group", "item", "a", "b")]
  for (m in c("haebara", "sl")) {
    expect_identical(link_groups(ip, method = m)$trend,
                     link_groups(ip_noc, method = m)$trend)
  }
})

test_that("item-varying guessing changes the solution under DIF", {
  ip <- ip_3pl(rep(c(0.05, 0.35), 12))
  ip$b[ip$group == 2][1:4] <- ip$b[ip$group == 2][1:4] + 0.5
  ip0 <- ip
  ip0$c <- 0
  for (m in c("haebara", "sl")) {
    l3 <- link_groups(ip, approach = "chain", method = m)
    l0 <- link_groups(ip0, approach = "chain", method = m)
    expect_true(all(l3$converged) && all(l0$converged))
    expect_gt(abs(l3$trend$mu[2] - l0$trend$mu[2]), 1e-6)
  }
})

test_that("pairwise joint criteria support the 3PL", {
  ip <- ip_3pl(0.2)
  lk <- link_groups(ip, approach = "joint", method = "haebara",
                    variant = "pairwise")
  expect_equal(lk$trend$mu[2], 0.3, tolerance = 1e-3)
  expect_equal(lk$trend$sigma[2], 1.1, tolerance = 1e-3)
  lk <- link_groups(ip, approach = "joint", method = "sl")
  expect_equal(lk$trend$mu[2], 0.3, tolerance = 1e-3)
})

test_that("guessing guards fire where closed forms are 2PL-based", {
  ip <- ip_3pl(0.2)
  expect_error(link_groups(ip, approach = "joint", method = "haebara"),
               "pairwise")
  for (m in c("haebara", "sl")) {
    lk <- link_groups(ip, approach = "chain", method = m)
    expect_error(linking_error(lk, method = "ajk"), "guessing")
    le <- linking_error(lk, method = "jackknife")
    expect_true(all(is.finite(le$le$le_mu[le$le$group == 2])))
  }
})

test_that("as_ipars validates the guessing column", {
  ip <- ip_3pl(1.2)
  expect_error(link_groups(ip), "\\[0, 1\\)")
})
