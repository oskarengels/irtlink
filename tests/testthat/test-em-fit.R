em_fit_sim <- function(seed = 31, N = 1500L, J = 12L, model = "2PL") {
  set.seed(seed)
  a <- if (model == "2PL") runif(J, 0.6, 1.8) else rep(1, J)
  b <- rnorm(J, 0, 1)
  theta <- rnorm(N)
  P <- plogis(sweep(outer(theta, a), 2, a * b, "-"))
  X <- matrix(rbinom(N * J, 1L, P), N, J)
  colnames(X) <- paste0("i", seq_len(J))
  list(X = X, a = a, b = b)
}

test_that("em_fit converges with a monotone log-likelihood (2PL)", {
  skip_on_cran()
  s <- em_fit_sim(31)
  fit <- em_fit(s$X, "2PL")
  expect_s3_class(fit, "irtlink_em")
  expect_true(fit$converged)
  expect_lt(fit$iter, fit$control$maxit)
  expect_true(all(diff(fit$loglik_trace) > -1e-6))
  expect_equal(fit$deviance, -2 * fit$loglik)
  expect_equal(fit$npar, 2L * ncol(s$X))
  expect_equal(fit$b, fit$nu / fit$a)
  expect_identical(fit$dat, {
    X <- s$X; storage.mode(X) <- "integer"; X
  })
})

test_that("em_fit recovers generating 2PL parameters approximately", {
  skip_on_cran()
  s <- em_fit_sim(32, N = 3000L)
  fit <- em_fit(s$X, "2PL")
  expect_gt(cor(fit$a, s$a), 0.9)
  expect_lt(mean(abs(fit$b - s$b)), 0.1)
})

test_that("em_fit 1PL fixes all discriminations at 1", {
  skip_on_cran()
  s <- em_fit_sim(33, model = "1PL")
  fit <- em_fit(s$X, "1PL")
  expect_true(all(fit$a == 1))
  expect_equal(fit$npar, ncol(s$X))
  expect_true(fit$converged)
})

test_that("pattern collapsing does not change the fit", {
  skip_on_cran()
  s <- em_fit_sim(34, N = 800L, J = 5L)  # <= 32 patterns -> collapsing on
  f1 <- em_fit(s$X, "2PL", control = list(collapse = TRUE))
  f0 <- em_fit(s$X, "2PL", control = list(collapse = FALSE))
  expect_equal(f1$a, f0$a, tolerance = 1e-8)
  expect_equal(f1$nu, f0$nu, tolerance = 1e-8)
  expect_equal(f1$loglik, f0$loglik, tolerance = 1e-8)
})

test_that("em_fit recovers a shifted distribution with fixed items", {
  skip_on_cran()
  set.seed(36)   # pure Hanson-1996 case: items known, distribution free
  J <- 12L; N <- 3000L
  a <- runif(J, 0.6, 1.8); b <- rnorm(J)
  th <- rnorm(N, 0.5, 1.2)
  P <- plogis(sweep(outer(th, a), 2, a * b, "-"))
  X <- matrix(rbinom(N * J, 1L, P), N, J)
  colnames(X) <- paste0("i", seq_len(J))
  bank <- data.frame(item = colnames(X), a = a, b = b,
                     stringsAsFactors = FALSE)
  fit <- em_fit(X, "2PL", fixed = bank, est_dist = TRUE)
  expect_true(fit$converged)
  expect_equal(unname(fit$dist[["mu"]]), 0.5, tolerance = 0.1)
  expect_equal(unname(fit$dist[["sigma"]]), 1.2, tolerance = 0.1)
  expect_equal(fit$a, a)                # fixed values untouched
  expect_equal(fit$b, b, tolerance = 1e-12)
  expect_equal(fit$npar, 2L)            # only mu, sigma estimated
})

test_that("em_fit estimates free items alongside fixed ones", {
  skip_on_cran()
  set.seed(37)
  J <- 10L; N <- 2000L
  a <- runif(J, 0.7, 1.6); b <- rnorm(J)
  th <- rnorm(N)
  P <- plogis(sweep(outer(th, a), 2, a * b, "-"))
  X <- matrix(rbinom(N * J, 1L, P), N, J)
  colnames(X) <- paste0("i", seq_len(J))
  bank <- data.frame(item = colnames(X)[1:6], a = a[1:6], b = b[1:6],
                     stringsAsFactors = FALSE)
  fit <- em_fit(X, "2PL", fixed = bank)
  expect_equal(fit$a[1:6], a[1:6])
  expect_equal(fit$b[1:6], b[1:6], tolerance = 1e-12)
  expect_gt(cor(fit$a[7:10], a[7:10]), 0.7)
  expect_lt(mean(abs(fit$b[7:10] - b[7:10])), 0.15)
  expect_equal(fit$npar, 2L * 4L)
  expect_equal(unname(fit$dist), c(0, 1))
})

test_that("em_fit estimates a 3PL with prior-stabilized guessing", {
  skip_on_cran()
  set.seed(38)
  J <- 15L; N <- 4000L
  a <- runif(J, 0.8, 1.8); b <- rnorm(J); cc <- runif(J, 0.05, 0.25)
  th <- rnorm(N)
  psi <- plogis(sweep(outer(th, a), 2, a * b, "-"))
  P <- sweep(sweep(psi, 2, 1 - cc, "*"), 2, cc, "+")
  X <- matrix(rbinom(N * J, 1L, P), N, J)
  colnames(X) <- paste0("i", seq_len(J))
  fit <- em_fit(X, "3PL")
  expect_true(fit$converged)
  expect_true(all(diff(fit$loglik_trace) > -1e-6))
  expect_length(fit$c, J)
  expect_true(all(fit$c >= 0.001 & fit$c <= 0.5))
  expect_gt(cor(fit$a, a), 0.8)
  expect_gt(cor(fit$b, b), 0.95)
  expect_lt(mean(abs(fit$c - cc)), 0.08)
  expect_equal(fit$npar, 3L * J)
})

test_that("3PL guards reject unsupported paths", {
  skip_on_cran()
  set.seed(39)
  X <- matrix(rbinom(600, 1L, 0.6), 60, 10,
              dimnames = list(NULL, paste0("i", 1:10)))
  bank <- data.frame(item = "i1", a = 1, b = 0)
  expect_error(em_fit(X, "3PL", fixed = bank), "3PL")
  expect_error(em_fit(X, "3PL", est_dist = TRUE), "3PL")
})

test_that("em_fit handles missing responses and unknown control errors", {
  skip_on_cran()
  s <- em_fit_sim(35)
  X <- s$X
  X[sample(length(X), floor(0.15 * length(X)))] <- NA
  fit <- em_fit(X, "2PL")
  expect_true(fit$converged)
  expect_error(em_fit(s$X, "2PL", control = list(bogus = 1)),
               "Unknown control option")
})
