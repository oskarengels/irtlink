test_that("the analytic score matches a numerical loglik gradient", {
  set.seed(91)
  d <- sim_trend_data(n_groups = 2, N = 400, I = 8, dif = "none")
  fit <- em_fit(as.matrix(as.data.frame(d[[1]])), "2PL")
  # loglik as a function of (a, nu), evaluated away from the optimum
  ll_at <- function(a, nu) {
    P <- em_irf_matrix(a, nu, fit$theta)
    ed <- em_estep_data(fit$dat, rep(1, nrow(fit$dat)))
    em_estep(ed, log(P), log(1 - P), log(fit$pi))$loglik
  }
  a0 <- fit$a * 0.9; nu0 <- fit$nu + 0.1
  sc <- em_score_ipars(fit, a0, nu0)
  eps <- 1e-5
  for (j in c(1L, 4L)) {
    ap <- a0; ap[j] <- ap[j] + eps
    am <- a0; am[j] <- am[j] - eps
    expect_equal(sc[2 * j - 1], (ll_at(ap, nu0) - ll_at(am, nu0)) / (2 * eps),
                 tolerance = 1e-4)
    vp <- nu0; vp[j] <- vp[j] + eps
    vm <- nu0; vm[j] <- vm[j] - eps
    expect_equal(sc[2 * j], (ll_at(a0, vp) - ll_at(a0, vm)) / (2 * eps),
                 tolerance = 1e-4)
  }
})

test_that("em_vcov_ipars returns a labeled a,b covariance", {
  set.seed(92)
  d <- sim_trend_data(n_groups = 2, N = 800, I = 8, dif = "none")
  X <- as.matrix(as.data.frame(d[[1]]))
  fit <- em_fit(X, "2PL")
  V <- em_vcov_ipars(fit)
  expect_equal(nrow(V), 2L * ncol(X))
  expect_true(all(grepl("_(a|b)$", rownames(V))))
  expect_true(all(diag(V) > 0))
  # SEs shrink roughly like 1/sqrt(N)
  fit2 <- em_fit(rbind(X, X, X, X), "2PL")
  V2 <- em_vcov_ipars(fit2)
  expect_equal(mean(sqrt(diag(V2)) / sqrt(diag(V))), 0.5, tolerance = 0.15)
})

test_that("1PL vcov has zero a rows and positive b variances", {
  set.seed(93)
  d <- sim_trend_data(n_groups = 2, N = 600, I = 8, dif = "none",
                      model = "1PL")
  fit <- em_fit(as.matrix(as.data.frame(d[[1]])), "1PL")
  V <- em_vcov_ipars(fit)
  a_rows <- grep("_a$", rownames(V))
  expect_true(all(V[a_rows, ] == 0))
  expect_true(all(diag(V)[grep("_b$", rownames(V))] > 0))
})

test_that("calibrate keep_vcov works for the em engine", {
  set.seed(94)
  d <- sim_trend_data(n_groups = 2, N = 500, I = 15, dif = "none")
  cal <- calibrate(d, engine = "em", keep_vcov = TRUE)
  expect_length(cal$vcov_ipars, 2L)
  expect_true(all(vapply(cal$vcov_ipars, is.matrix, logical(1))))
})

test_that("the 3PL score matches a numerical penalized-loglik gradient", {
  set.seed(95)
  d <- sim_trend_data(n_groups = 2, N = 500, I = 6, dif = "none")
  fit <- em_fit(as.matrix(as.data.frame(d[[1]])), "3PL")
  # Bayes-modal objective: marginal loglik plus the beta log-prior on c,
  # evaluated away from the optimum so no score component is trivially 0.
  pl_at <- function(a, nu, gu) {
    P <- em_irf_matrix(a, nu, fit$theta, guess = gu)
    ed <- em_estep_data(fit$dat, rep(1, nrow(fit$dat)))
    ll <- em_estep(ed, log(P), log(1 - P), log(fit$pi))$loglik
    pc <- fit$control$prior_c
    ll + sum((pc[1] - 1) * log(gu) + (pc[2] - 1) * log(1 - gu))
  }
  a0 <- fit$a * 0.95; nu0 <- fit$nu + 0.05
  g0 <- pmin(pmax(fit$c, 0.05), 0.3)
  sc <- em_score_ipars(fit, a0, nu0, gu = g0)
  eps <- 1e-5
  for (j in c(1L, 3L)) {
    ap <- a0; ap[j] <- ap[j] + eps
    am <- a0; am[j] <- am[j] - eps
    expect_equal(sc[3 * j - 2],
                 (pl_at(ap, nu0, g0) - pl_at(am, nu0, g0)) / (2 * eps),
                 tolerance = 1e-4)
    vp <- nu0; vp[j] <- vp[j] + eps
    vm <- nu0; vm[j] <- vm[j] - eps
    expect_equal(sc[3 * j - 1],
                 (pl_at(a0, vp, g0) - pl_at(a0, vm, g0)) / (2 * eps),
                 tolerance = 1e-4)
    gp <- g0; gp[j] <- gp[j] + eps
    gm <- g0; gm[j] <- gm[j] - eps
    expect_equal(sc[3 * j],
                 (pl_at(a0, nu0, gp) - pl_at(a0, nu0, gm)) / (2 * eps),
                 tolerance = 1e-4)
  }
})

test_that("em_vcov_ipars returns a labeled a,b,c covariance for the 3PL", {
  set.seed(96)
  d <- sim_trend_data(n_groups = 2, N = 1000, I = 8, dif = "none")
  X <- as.matrix(as.data.frame(d[[1]]))
  fit <- em_fit(X, "3PL")
  V <- em_vcov_ipars(fit)
  expect_equal(nrow(V), 3L * ncol(X))
  expect_identical(rownames(V)[1:3], paste0(fit$item[1], c("_a", "_b", "_c")))
  expect_lt(max(abs(V - t(V))), 1e-10)
  expect_true(all(diag(V) > 0))
  ev <- eigen(V, symmetric = TRUE, only.values = TRUE)$values
  expect_gt(min(ev), -1e-8 * max(ev))
})

test_that("calibrate keep_vcov works for the 3PL and feeds linking_error", {
  set.seed(97)
  d <- sim_trend_data(n_groups = 2, N = 800, I = 10, overlap = 1, dif = "none")
  cal <- calibrate(d, model = "3PL", keep_vcov = TRUE)
  expect_length(cal$vcov_ipars, 2L)
  expect_true(all(grepl("_(a|b|c)$", rownames(cal$vcov_ipars[[1]]))))
  # le_vcov_matrix pulls the named a,b block out of the 3 x 3 item blocks
  lk <- link_chain(cal, method = "mgm")
  le <- linking_error(lk, method = "jackknife_bc", vcov = cal$vcov_ipars)$le
  expect_true(all(is.finite(le$se_mu)))
  expect_gt(max(le$se_mu), 0)
})
