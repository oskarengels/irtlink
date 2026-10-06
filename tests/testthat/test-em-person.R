em_person_sim <- function(seed = 111, N = 500L, J = 12L) {
  set.seed(seed)
  a <- runif(J, 0.7, 1.6); b <- rnorm(J)
  th <- rnorm(N)
  P <- plogis(sweep(outer(th, a), 2, a * b, "-"))
  X <- matrix(rbinom(N * J, 1L, P), N, J)
  colnames(X) <- paste0("i", seq_len(J))
  storage.mode(X) <- "integer"
  list(X = X, a = a, nu = a * b, th = th)
}

test_that("em_eap_core matches a direct posterior computation", {
  s <- em_person_sim()
  theta <- em_theta_grid(41L, c(-5, 5))
  pi_k <- dnorm_discrete(theta, 0, 1)
  eap <- em_eap_core(s$X, s$a, s$nu, theta, pi_k)
  expect_equal(nrow(eap), nrow(s$X))
  for (i in c(1L, 17L, 300L)) {
    lk <- log(pi_k)
    for (j in seq_len(ncol(s$X))) {
      if (is.na(s$X[i, j])) next
      p <- plogis(s$a[j] * theta - s$nu[j])
      lk <- lk + if (s$X[i, j] == 1L) log(p) else log(1 - p)
    }
    w <- exp(lk - max(lk)); w <- w / sum(w)
    expect_equal(eap$est[i], sum(w * theta), tolerance = 1e-10)
    expect_equal(eap$se[i],
                 sqrt(sum(w * theta^2) - sum(w * theta)^2),
                 tolerance = 1e-10)
  }
  expect_gt(cor(eap$est, s$th), 0.8)
  expect_true(all(eap$se > 0))
})

test_that("em_eap_core handles missing responses per person", {
  s <- em_person_sim(112)
  X <- s$X
  X[1, ] <- NA_integer_        # no data: posterior = prior
  theta <- em_theta_grid(41L, c(-5, 5))
  pi_k <- dnorm_discrete(theta, 0, 1)
  eap <- em_eap_core(X, s$a, s$nu, theta, pi_k)
  expect_equal(eap$est[1], sum(pi_k * theta), tolerance = 1e-10)
})

test_that("em_wle_core solves the Warm score equation per person", {
  s <- em_person_sim(113)
  wle <- em_wle_core(s$X, s$a, s$nu)
  expect_equal(nrow(wle), nrow(s$X))
  for (i in c(2L, 50L, 400L)) {
    th <- wle$est[i]
    p <- plogis(s$a * th - s$nu)
    obs <- !is.na(s$X[i, ])
    info <- sum((s$a^2 * p * (1 - p))[obs])
    jterm <- sum((s$a^3 * p * (1 - p) * (1 - 2 * p))[obs])
    score <- sum((s$a * (s$X[i, ] - p))[obs]) + jterm / (2 * info)
    expect_lt(abs(score), 1e-6)
    expect_equal(wle$se[i], 1 / sqrt(info), tolerance = 1e-8)
  }
  expect_gt(cor(wle$est, s$th), 0.8)
})

test_that("em_wle_core stays finite for perfect and zero scores", {
  s <- em_person_sim(114, N = 4L, J = 10L)
  X <- s$X
  X[1, ] <- 1L
  X[2, ] <- 0L
  wle <- em_wle_core(X, s$a, s$nu)
  expect_true(all(is.finite(wle$est)))
  expect_true(all(is.finite(wle$se)))
  expect_gt(wle$est[1], 1)
  expect_lt(wle$est[2], -1)
})

test_that("person_scores works for separate em calibrations", {
  set.seed(115)
  d <- sim_trend_data(n_groups = 2, N = 300, I = 15, dif = "none")
  cal <- calibrate(d, engine = "em")
  for (m in c("eap", "wle")) {
    ps <- person_scores(cal, method = m)
    expect_named(ps, c("group", "person", "est", "se"))
    expect_equal(nrow(ps), sum(vapply(d, nrow, integer(1))))
    expect_true(all(is.finite(ps$est)))
  }
})

test_that("person_scores works for concurrent em calibrations", {
  set.seed(116)
  d <- sim_trend_data(n_groups = 2, N = 300, I = 15, dif = "none")
  cal <- calibrate(d, calibration = "concurrent", engine = "em")
  ps <- person_scores(cal, method = "eap")
  expect_equal(sort(unique(ps$group)), 1:2)
  # group-2 EAPs shrink toward the group-2 prior mean, not toward 0
  expect_equal(mean(ps$est[ps$group == 2]), cal$trend$mu[2],
               tolerance = 0.15)
})

test_that("person_scores rejects unsupported inputs", {
  set.seed(117)
  d <- sim_trend_data(n_groups = 2, N = 200, I = 15, dif = "none")
  cal_nm <- calibrate(d, keep_models = FALSE)
  expect_error(person_scores(cal_nm), "keep_models")
})

