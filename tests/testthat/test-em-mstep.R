em_mstep_case <- function(seed) {
  s <- em_ref_setup(seed = seed, N = 300L, J = 5L, K = 21L, na_frac = 0.05)
  es <- em_estep_ref(s$X, s$P, s$pi)
  list(s = s, es = es)
}

test_that("em_mstep_items_cpp matches the optim reference (2PL)", {
  cs <- em_mstep_case(21)
  J <- length(cs$s$a)
  a0 <- rep(1, J); nu0 <- rep(0, J)
  cpp <- em_mstep_items_cpp(cs$s$theta, cs$es$njk, cs$es$rjk, a0, nu0,
                            rep(TRUE, J), rep(0.01, J), rep(10, J),
                            rep(-10, J), rep(10, J), 100L, 1e-10)
  ref <- em_mstep_ref(cs$s$theta, cs$es$njk, cs$es$rjk, a0, nu0,
                      rep(TRUE, J))
  expect_equal(cpp$a, ref$a, tolerance = 1e-4)
  expect_equal(cpp$nu, ref$nu, tolerance = 1e-4)
})

test_that("the analytic score is zero at the kernel solution", {
  cs <- em_mstep_case(22)
  J <- length(cs$s$a)
  cpp <- em_mstep_items_cpp(cs$s$theta, cs$es$njk, cs$es$rjk,
                            rep(1, J), rep(0, J), rep(TRUE, J),
                            rep(0.01, J), rep(10, J), rep(-10, J),
                            rep(10, J), 100L, 1e-12)
  for (j in seq_len(J)) {
    p <- plogis(cpp$a[j] * cs$s$theta - cpp$nu[j])
    resid <- cs$es$rjk[j, ] - cs$es$njk[j, ] * p
    expect_lt(abs(sum(resid * cs$s$theta)), 1e-6)
    expect_lt(abs(sum(resid)), 1e-6)
  }
})

test_that("est_a = FALSE fixes a (1PL) and only updates nu", {
  cs <- em_mstep_case(23)
  J <- length(cs$s$a)
  cpp <- em_mstep_items_cpp(cs$s$theta, cs$es$njk, cs$es$rjk,
                            rep(1, J), rep(0, J), rep(FALSE, J),
                            rep(0.01, J), rep(10, J), rep(-10, J),
                            rep(10, J), 100L, 1e-10)
  expect_equal(cpp$a, rep(1, J))
  ref <- em_mstep_ref(cs$s$theta, cs$es$njk, cs$es$rjk,
                      rep(1, J), rep(0, J), rep(FALSE, J))
  expect_equal(cpp$nu, ref$nu, tolerance = 1e-4)
})

test_that("em_mstep_item_3pl recovers (a, nu, c) from noiseless counts", {
  set.seed(25)
  theta <- em_theta_grid(31L, c(-4, 4))
  njk <- 5000 * dnorm_discrete(theta, 0, 1)
  a_true <- 1.3; v_true <- 0.4; c_true <- 0.2
  P <- c_true + (1 - c_true) * plogis(a_true * theta - v_true)
  rjk <- njk * P
  res <- em_mstep_item_3pl(theta, njk, rjk, a0 = 1, nu0 = 0, c0 = 0.15,
                           prior_c = c(1, 1))   # flat prior: pure ML
  expect_equal(res$a, a_true, tolerance = 1e-3)
  expect_equal(res$nu, v_true, tolerance = 1e-3)
  expect_equal(res$c, c_true, tolerance = 1e-3)
})

test_that("the beta prior pulls c toward its mode on weak data", {
  set.seed(26)
  theta <- em_theta_grid(31L, c(-4, 4))
  njk <- 30 * dnorm_discrete(theta, 0, 1)     # very little data
  P <- 0.05 + 0.95 * plogis(1.2 * theta)
  rjk <- njk * P
  flat <- em_mstep_item_3pl(theta, njk, rjk, 1, 0, 0.15, prior_c = c(1, 1))
  informative <- em_mstep_item_3pl(theta, njk, rjk, 1, 0, 0.15,
                                   prior_c = c(5, 17))
  expect_gt(informative$c, flat$c)   # pulled toward the 0.2 prior mode
})

test_that("the 3PL analytic gradient matches finite differences", {
  set.seed(27)
  theta <- em_theta_grid(21L, c(-4, 4))
  njk <- 200 * dnorm_discrete(theta, 0.2, 1.1)
  rjk <- njk * (0.15 + 0.85 * plogis(1.1 * theta - 0.3))
  par <- c(0.9, 0.2, 0.12); prior_c <- c(5, 17)
  gr <- em_mstep_item_3pl_gradient(par, theta, njk, rjk, prior_c)
  f <- function(p) em_mstep_item_3pl_objective(p, theta, njk, rjk, prior_c)
  eps <- 1e-6
  for (j in 1:3) {
    pp <- par; pp[j] <- pp[j] + eps
    pm <- par; pm[j] <- pm[j] - eps
    expect_equal(gr[j], (f(pp) - f(pm)) / (2 * eps), tolerance = 1e-4)
  }
})

test_that("degenerate items stay within bounds", {
  theta <- em_theta_grid(11L, c(-4, 4))
  njk <- matrix(10, 1L, 11L)
  rjk <- njk  # all responses correct -> nu runs to the lower bound
  cpp <- em_mstep_items_cpp(theta, njk, rjk, 1, 0, TRUE,
                            0.01, 10, -10, 10, 200L, 1e-10)
  expect_gte(cpp$a, 0.01); expect_lte(cpp$a, 10)
  expect_gte(cpp$nu, -10); expect_lte(cpp$nu, 10)
})
