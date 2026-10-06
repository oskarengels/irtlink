# Pure-R reference implementations of the EM kernels. Used only to
# verify the C++ kernels; intentionally simple and slow.

em_estep_ref <- function(X, P, pi, weights = rep(1, nrow(X))) {
  N <- nrow(X); J <- ncol(X); K <- length(pi)
  njk <- matrix(0, J, K); rjk <- matrix(0, J, K); ll <- 0
  for (i in seq_len(N)) {
    lk <- log(pi)
    for (j in seq_len(J)) {
      y <- X[i, j]
      if (is.na(y)) next
      lk <- lk + if (y == 1L) log(P[j, ]) else log(1 - P[j, ])
    }
    m <- max(lk); ew <- exp(lk - m); s <- sum(ew)
    ll <- ll + weights[i] * (m + log(s))
    w <- weights[i] * ew / s
    for (j in seq_len(J)) {
      y <- X[i, j]
      if (is.na(y)) next
      njk[j, ] <- njk[j, ] + w
      if (y == 1L) rjk[j, ] <- rjk[j, ] + w
    }
  }
  list(njk = njk, rjk = rjk, loglik = ll)
}

em_mstep_ref <- function(theta, njk, rjk, a, nu, est_a,
                         lower = c(0.01, -10), upper = c(10, 10)) {
  for (j in seq_along(a)) {
    nll <- function(par) {
      aj <- if (est_a[j]) par[1] else a[j]
      vj <- if (est_a[j]) par[2] else par[1]
      p <- plogis(aj * theta - vj)
      p <- pmin(pmax(p, 1e-10), 1 - 1e-10)
      -sum(rjk[j, ] * log(p) + (njk[j, ] - rjk[j, ]) * log(1 - p))
    }
    if (est_a[j]) {
      o <- stats::optim(c(a[j], nu[j]), nll, method = "L-BFGS-B",
                        lower = lower, upper = upper)
      a[j] <- o$par[1]; nu[j] <- o$par[2]
    } else {
      o <- stats::optim(nu[j], nll, method = "L-BFGS-B",
                        lower = lower[2], upper = upper[2])
      nu[j] <- o$par[1]
    }
  }
  list(a = a, nu = nu)
}

em_ref_setup <- function(seed = 1, N = 40L, J = 6L, K = 21L,
                         na_frac = 0.1) {
  set.seed(seed)
  a <- runif(J, 0.7, 1.8)
  nu <- rnorm(J, 0, 1)
  theta_true <- rnorm(N)
  P_true <- plogis(outer(a, theta_true) - nu)
  X <- matrix(rbinom(N * J, 1L, t(P_true)), nrow = N, ncol = J)
  if (na_frac > 0) X[sample(length(X), floor(na_frac * length(X)))] <- NA
  storage.mode(X) <- "integer"
  theta <- em_theta_grid(K, c(-4, 4))
  pi <- dnorm_discrete(theta, 0, 1)
  P <- em_irf_matrix(a, nu, theta)
  list(X = X, theta = theta, pi = pi, P = P, a = a, nu = nu)
}
