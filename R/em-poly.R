# R/em-poly.R
# generalized partial credit model (GPCM, Muraki, 1992) and partial
# credit model (PCM, Masters, 1982) for the internal EM engine

# category probabilities of one item on a theta grid, d holds the
# cumulative intercepts, eta_m = m a theta - d_m with d_0 = 0
em_poly_probs <- function(a, d, theta) {
  K <- length(d)
  eta <- outer(theta, 0:K) * a
  eta <- sweep(eta, 2L, c(0, d), "-")
  eta <- eta - apply(eta, 1L, max)
  P <- exp(eta)
  P <- P / rowSums(P)
  pmin(pmax(P, 1e-10), 1)
}

# log category probabilities per item and node, zeros where the
# category does not exist for the item
em_poly_logP <- function(a, dmat, K, theta) {
  J <- length(a)
  Pl <- lapply(seq_len(J), function(j)
    em_poly_probs(a[j], dmat[j, seq_len(K[j])], theta))
  lapply(0:max(K), function(m) {
    L <- matrix(0, J, length(theta))
    for (j in seq_len(J)) if (m <= K[j]) L[j, ] <- log(Pl[[j]][, m + 1])
    L
  })
}

# response masks per category for the polytomous E step
em_estep_data_poly <- function(X, w, K) {
  Kmax <- max(K)
  Obs <- (!is.na(X)) * 1
  Ym <- lapply(0:Kmax, function(m) {
    Y <- (!is.na(X) & X == m) * 1
    storage.mode(Y) <- "double"
    Y
  })
  list(Ym = Ym, Obs = Obs, w = w, N = nrow(X), Kmax = Kmax)
}

# polytomous E step, logP[[m + 1]] holds log P(X = m) per item and
# node with zeros where the category does not exist for the item
em_estep_poly <- function(ed, logP, log_pi, return_posterior = FALSE) {
  LL <- ed$Ym[[1]] %*% logP[[1]]
  for (m in seq_len(ed$Kmax)) LL <- LL + ed$Ym[[m + 1]] %*% logP[[m + 1]]
  LL <- sweep(LL, 2L, log_pi, "+")
  mx <- LL[cbind(seq_len(ed$N), max.col(LL, ties.method = "first"))]
  W <- exp(LL - mx)
  s <- rowSums(W)
  loglik <- sum(ed$w * (mx + log(s)))
  W <- (ed$w / s) * W
  njk <- crossprod(ed$Obs, W)
  rjkm <- lapply(ed$Ym, function(Y) crossprod(Y, W))
  out <- list(njk = njk, rjkm = rjkm, nk = colSums(W), loglik = loglik)
  if (return_posterior) out$posterior <- W / ed$w
  out
}

# item M step, Newton steps with step halving on the concave
# complete-data loglik in (a, d_1, ..., d_K)
em_mstep_item_poly <- function(theta, njk, rjm, a0, d0, est_a,
                               lower_a = 0.1, upper_a = 10,
                               maxit = 5L, tol = 1e-8) {
  K <- length(d0)
  mvec <- 0:K
  phi <- c(a0, d0)
  cll <- function(phi) {
    P <- em_poly_probs(phi[1], phi[-1], theta)
    sum(rjm * log(P))
  }
  f0 <- cll(phi)
  for (it in seq_len(maxit)) {
    a <- phi[1]
    P <- em_poly_probs(a, phi[-1], theta)
    Ebar <- as.numeric(P %*% mvec)
    Vm <- as.numeric(P %*% mvec^2) - Ebar^2
    rb <- as.numeric(rjm %*% mvec)
    g_a <- sum(theta * (rb - njk * Ebar))
    g_d <- colSums(njk * P[, -1, drop = FALSE]) -
      colSums(rjm[, -1, drop = FALSE])
    H <- matrix(0, K + 1, K + 1)
    H[1, 1] <- -sum(theta^2 * njk * Vm)
    for (m in seq_len(K)) {
      H[1, m + 1] <- H[m + 1, 1] <-
        sum(theta * njk * P[, m + 1] * (m - Ebar))
      for (l in m:K) {
        v <- sum(njk * P[, m + 1] * (P[, l + 1] - (m == l)))
        H[m + 1, l + 1] <- H[l + 1, m + 1] <- v
      }
    }
    g <- c(g_a, g_d)
    if (est_a) {
      delta <- tryCatch(solve(-H, g), error = function(e) NULL)
    } else {
      delta <- tryCatch(c(0, solve(-H[-1, -1, drop = FALSE], g_d)),
                        error = function(e) NULL)
    }
    if (is.null(delta)) break
    stp <- 1
    repeat {
      cand <- phi + stp * delta
      cand[1] <- min(max(cand[1], lower_a), upper_a)
      f1 <- cll(cand)
      if (is.finite(f1) && f1 >= f0 - 1e-12) break
      stp <- stp / 2
      if (stp < 1e-4) {
        cand <- phi
        f1 <- f0
        break
      }
    }
    moved <- max(abs(cand - phi))
    phi <- cand
    f0 <- f1
    if (moved < tol) break
  }
  list(a = phi[1], d = phi[-1])
}

# EM loop for one group under a fixed N(0,1) distribution
em_fit_poly <- function(X, model = c("GPCM", "PCM"), control = list(),
                        pweights = NULL) {
  model <- match.arg(model)
  ctrl <- em_control(control)
  X <- as.matrix(X)
  storage.mode(X) <- "integer"
  if (is.null(colnames(X))) colnames(X) <- paste0("I", seq_len(ncol(X)))
  no_obs <- colSums(!is.na(X)) == 0L
  if (any(no_obs)) {
    stop("Item(s) ", paste(colnames(X)[no_obs], collapse = ", "),
         " have no observed responses in this group.", call. = FALSE)
  }
  J <- ncol(X)
  K <- apply(X, 2L, max, na.rm = TRUE)
  if (any(K < 1)) {
    stop("Item(s) ", paste(colnames(X)[K < 1], collapse = ", "),
         " have a single observed category.", call. = FALSE)
  }
  est_a <- model == "GPCM"
  lower_a <- if (!is.null(ctrl$lower_a)) ctrl$lower_a else 0.1
  theta <- em_theta_grid(ctrl$theta_nodes, ctrl$theta_range)
  # starts, a = 1 and equal thresholds from the scaled item means
  pbar <- vapply(seq_len(J), function(j)
    mean(X[, j], na.rm = TRUE) / K[j], numeric(1))
  pbar <- pmin(pmax(pbar, 0.05), 0.95)
  b_start <- -stats::qlogis(pbar)
  a <- rep(1, J)
  dmat <- matrix(NA_real_, J, max(K))
  for (j in seq_len(J)) dmat[j, seq_len(K[j])] <- cumsum(rep(b_start[j], K[j]))
  cp <- if (ctrl$collapse) em_collapse_patterns(X, pweights)
        else list(X = X, w = if (is.null(pweights)) rep(1, nrow(X))
                             else pweights)
  ed <- em_estep_data_poly(cp$X, cp$w, K)
  pi_k <- dnorm_discrete(theta, mean = 0, sd = 1)
  log_pi <- log(pi_k)
  ll_trace <- numeric(0)
  converged <- FALSE
  iter <- 0L
  for (s in seq_len(ctrl$maxit)) {
    es <- em_estep_poly(ed, em_poly_logP(a, dmat, K, theta), log_pi)
    if (length(ll_trace) > 0 &&
        es$loglik < ll_trace[length(ll_trace)] -
          1e-8 * (1 + abs(es$loglik))) {
      warning("EM log-likelihood decreased at iteration ", s, ".",
              call. = FALSE)
    }
    ll_trace <- c(ll_trace, es$loglik)
    old <- c(a, dmat[!is.na(dmat)])
    for (j in seq_len(J)) {
      rjm <- vapply(0:K[j], function(m) es$rjkm[[m + 1]][j, ],
                    numeric(length(theta)))
      up <- em_mstep_item_poly(theta, es$njk[j, ], rjm, a[j],
                               dmat[j, seq_len(K[j])], est_a,
                               lower_a = lower_a,
                               maxit = ctrl$inner_maxit,
                               tol = ctrl$inner_tol)
      a[j] <- up$a
      dmat[j, seq_len(K[j])] <- up$d
    }
    delta <- max(abs(c(a, dmat[!is.na(dmat)]) - old))
    iter <- s
    if (ctrl$verbose) {
      cat(sprintf("EM iter %d: loglik %.6f, max change %.2e\n",
                  s, es$loglik, delta))
    }
    if (delta < ctrl$conv) {
      converged <- TRUE
      break
    }
  }
  es <- em_estep_poly(ed, em_poly_logP(a, dmat, K, theta), log_pi)
  # thresholds b_v = (d_v - d_{v-1}) / a, location b = mean(b_v)
  thr <- matrix(NA_real_, J, max(K))
  for (j in seq_len(J)) {
    dj <- dmat[j, seq_len(K[j])]
    thr[j, seq_len(K[j])] <- diff(c(0, dj)) / a[j]
  }
  b <- rowMeans(thr, na.rm = TRUE)
  tau <- thr - b
  structure(
    list(item = colnames(X), a = as.numeric(a), b = as.numeric(b),
         thresholds = thr, tau = tau, d = dmat, K = as.integer(K),
         model = model,
         theta = theta, pi = pi_k,
         dist = c(mu = 0, sigma = 1),
         loglik = es$loglik,
         loglik_trace = ll_trace, deviance = -2 * es$loglik,
         npar = sum(K) + est_a * J,
         dat = X, pweights = pweights,
         converged = converged, iter = iter, control = ctrl),
    class = "irtlink_em"
  )
}

# EAP and posterior SD per person for a polytomous fit
em_eap_core_poly <- function(X, a, dmat, K, theta, pi) {
  X <- as.matrix(X)
  storage.mode(X) <- "integer"
  ed <- em_estep_data_poly(X, rep(1, nrow(X)), K)
  es <- em_estep_poly(ed, em_poly_logP(a, dmat, K, theta), log(pi),
                      return_posterior = TRUE)
  W <- es$posterior
  est <- as.numeric(W %*% theta)
  m2 <- as.numeric(W %*% theta^2)
  data.frame(est = est, se = sqrt(pmax(m2 - est^2, 0)))
}
