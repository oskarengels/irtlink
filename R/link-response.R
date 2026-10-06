# R/link-response.R
# two-group response function fit for one chain step
# (Haebara, 1980. Stocking & Lord, 1983)

loss_rho <- function(e, pow, eps) {
  if (pow == 2) return(e^2)
  if (pow == 0) return(log(e^2 + eps))
  (e^2 + eps)^(pow / 2)
}

loss_rho1 <- function(e, pow, eps) {
  if (pow == 2) return(2 * e)
  if (pow == 0) return(2 * e / (e^2 + eps))
  pow * e * (e^2 + eps)^(pow / 2 - 1)
}

# value and gradient in (mu, log sigma) as closures for optim(),
# c0/c1 are the 3PL guessing parameters, invariant under the
# transformation (Kolen & Brennan, 2014)
fit_link_response <- function(a0, b0, a1, b1, c0 = NULL, c1 = NULL,
                           method = c("Hae", "SL"),
                           type = c("asymm", "symm"), pow = 2,
                           eps = 1e-3, theta, wgt, maxit = 1000L) {
  method <- match.arg(method)
  type <- match.arg(type)
  if (is.null(c0)) c0 <- rep(0, length(a0))
  if (is.null(c1)) c1 <- rep(0, length(a1))
  with_guess <- function(L, cc) {
    sweep(sweep(L, 2L, 1 - cc, "*"), 2L, cc, "+")
  }
  P2 <- with_guess(vapply(seq_along(a1), function(i)
    stats::plogis(a1[i] * (theta - b1[i])), numeric(length(theta))), c1)
  P1f <- with_guess(vapply(seq_along(a0), function(i)
    stats::plogis(a0[i] * (theta - b0[i])), numeric(length(theta))), c0)

  pieces <- function(mu, sigma) {
    th1 <- sigma * theta + mu
    L1 <- vapply(seq_along(a0), function(i)
      stats::plogis(a0[i] * (th1 - b0[i])), numeric(length(theta)))
    d1 <- with_guess(L1, c0) - P2
    u1 <- sweep(L1 * (1 - L1), 2L, (1 - c0) * a0, "*")
    out <- list(d1 = d1, du1_mu = u1, du1_si = u1 * theta)
    if (type == "symm") {
      psi2 <- (theta - mu) / sigma
      L2 <- vapply(seq_along(a1), function(i)
        stats::plogis(a1[i] * (psi2 - b1[i])), numeric(length(theta)))
      w2 <- sweep(L2 * (1 - L2), 2L, (1 - c1) * a1, "*")
      out$d2 <- P1f - with_guess(L2, c1)
      out$du2_mu <- w2 / sigma
      out$du2_si <- w2 * (theta - mu) / sigma^2
    }
    out
  }

  obj <- function(par) {
    mu <- par[1L]
    sigma <- exp(par[2L])
    pc <- pieces(mu, sigma)
    if (method == "Hae") {
      val <- sum(wgt * rowSums(loss_rho(pc$d1, pow, eps)))
      if (type == "symm") {
        val <- val + sum(wgt * rowSums(loss_rho(pc$d2, pow, eps)))
      }
    } else {
      val <- sum(wgt * loss_rho(rowSums(pc$d1), pow, eps))
      if (type == "symm") {
        val <- val + sum(wgt * loss_rho(rowSums(pc$d2), pow, eps))
      }
    }
    val
  }

  grd <- function(par) {
    mu <- par[1L]
    sigma <- exp(par[2L])
    pc <- pieces(mu, sigma)
    if (method == "Hae") {
      r1 <- loss_rho1(pc$d1, pow, eps)
      g_mu <- sum(wgt * rowSums(r1 * pc$du1_mu))
      g_si <- sum(wgt * rowSums(r1 * pc$du1_si))
      if (type == "symm") {
        r2 <- loss_rho1(pc$d2, pow, eps)
        g_mu <- g_mu + sum(wgt * rowSums(r2 * pc$du2_mu))
        g_si <- g_si + sum(wgt * rowSums(r2 * pc$du2_si))
      }
    } else {
      r1 <- loss_rho1(rowSums(pc$d1), pow, eps)
      g_mu <- sum(wgt * r1 * rowSums(pc$du1_mu))
      g_si <- sum(wgt * r1 * rowSums(pc$du1_si))
      if (type == "symm") {
        r2 <- loss_rho1(rowSums(pc$d2), pow, eps)
        g_mu <- g_mu + sum(wgt * r2 * rowSums(pc$du2_mu))
        g_si <- g_si + sum(wgt * r2 * rowSums(pc$du2_si))
      }
    }
    c(g_mu, g_si * sigma)
  }

  sigma0 <- exp(mean(log(a1)) - mean(log(a0)))
  mu0 <- mean(b0) - sigma0 * mean(b1)
  opt <- stats::optim(c(mu0, log(sigma0)), fn = obj, gr = grd,
                      method = "BFGS",
                      control = list(maxit = maxit, reltol = 1e-14))
  list(mu = opt$par[1L], sigma = exp(opt$par[2L]),
       converged = opt$convergence == 0, value = opt$value)
}

# polytomous criteria for one chain step, Haebara over the category
# response functions and Stocking-Lord over the expected score
# functions (Kim & Lee, 2006)
fit_link_response_poly <- function(a0, b0, tau0, a1, b1, tau1,
                                   method = c("Hae", "SL"),
                                   type = c("asymm", "symm"), pow = 2,
                                   eps = 1e-3, theta, wgt,
                                   maxit = 1000L) {
  method <- match.arg(method)
  type <- match.arg(type)
  K0 <- rowSums(!is.na(tau0))
  K1 <- rowSums(!is.na(tau1))
  if (!identical(K0, K1)) {
    bad <- which(K0 != K1)
    stop("Common item(s) ", paste(bad, collapse = ", "), " have a ",
         "different number of observed categories in the two groups.",
         call. = FALSE)
  }
  thr0 <- lapply(seq_along(a0), function(j)
    b0[j] + tau0[j, seq_len(K0[j])])
  thr1 <- lapply(seq_along(a1), function(j)
    b1[j] + tau1[j, seq_len(K1[j])])
  cat_probs <- function(a_j, thr_j, th) {
    em_poly_probs(a_j, a_j * cumsum(thr_j), th)
  }
  crit_side <- function(th_from, th_to) {
    if (method == "Hae") {
      val <- 0
      for (j in seq_along(a0)) {
        e <- cat_probs(a0[j], thr0[[j]], th_from)[, -1, drop = FALSE] -
          cat_probs(a1[j], thr1[[j]], th_to)[, -1, drop = FALSE]
        val <- val + sum(wgt * rowSums(loss_rho(e, pow, eps)))
      }
      return(val)
    }
    T0 <- 0
    T1 <- 0
    for (j in seq_along(a0)) {
      T0 <- T0 + as.numeric(cat_probs(a0[j], thr0[[j]], th_from) %*%
                              (0:K0[j]))
      T1 <- T1 + as.numeric(cat_probs(a1[j], thr1[[j]], th_to) %*%
                              (0:K1[j]))
    }
    sum(wgt * loss_rho(T0 - T1, pow, eps))
  }
  obj <- function(par) {
    mu <- par[1L]
    sigma <- exp(par[2L])
    val <- crit_side(sigma * theta + mu, theta)
    if (type == "symm") {
      val <- val + crit_side(theta, (theta - mu) / sigma)
    }
    val
  }
  sigma0 <- exp(mean(log(a1)) - mean(log(a0)))
  mu0 <- mean(b0) - sigma0 * mean(b1)
  opt <- stats::optim(c(mu0, log(sigma0)), fn = obj, method = "BFGS",
                      control = list(maxit = maxit, reltol = 1e-12))
  list(mu = opt$par[1L], sigma = exp(opt$par[2L]),
       converged = opt$convergence == 0, value = opt$value)
}
