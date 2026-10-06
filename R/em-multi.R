# R/em-multi.R
# multiple group estimation for the internal EM engine

# prepare a data_list for the multiple group EM
em_multi_data <- function(data_list, collapse = TRUE, pweights = NULL) {
  Xs <- lapply(data_list, function(d) {
    X <- as.matrix(d)
    storage.mode(X) <- "integer"
    X
  })
  if (any(vapply(Xs, function(X) is.null(colnames(X)), logical(1))))
    stop("All groups need item column names.", call. = FALSE)
  pool <- sort(unique(unlist(lapply(Xs, colnames))))
  groups <- lapply(seq_along(Xs), function(t) {
    X <- Xs[[t]]
    wt <- if (is.null(pweights)) NULL else pweights[[t]]
    cp <- if (collapse) em_collapse_patterns(X, wt)
          else list(X = X, w = if (is.null(wt)) rep(1, nrow(X)) else wt)
    list(ed = em_estep_data(cp$X, cp$w),
         idx = match(colnames(X), pool),
         N = sum(cp$w))
  })
  list(groups = groups, pool = pool, n_groups = length(Xs), J = length(pool))
}

# per-group E step over the item union
em_estep_multi <- function(md, a, nu, g, theta, pis, h = NULL) {
  if (is.null(h)) h <- matrix(0, md$J, md$n_groups)
  K <- length(theta)
  groups <- vector("list", md$n_groups)
  loglik <- 0
  for (t in seq_len(md$n_groups)) {
    wv <- md$groups[[t]]
    idx <- wv$idx
    P <- em_irf_matrix(a[idx] + h[idx, t], nu[idx] + g[idx, t], theta)
    es <- em_estep(wv$ed, log(P), log(1 - P), log(pis[[t]]))
    njk <- matrix(0, md$J, K)
    rjk <- matrix(0, md$J, K)
    njk[idx, ] <- es$njk
    rjk[idx, ] <- es$rjk
    groups[[t]] <- list(njk = njk, rjk = rjk, nk = es$nk)
    loglik <- loglik + es$loglik
  }
  list(groups = groups, loglik = loglik)
}

# normal moment update on the node totals (Hanson, 1996)
em_update_normal <- function(nk, theta, N) {
  mu <- sum(theta * nk) / N
  sig2 <- sum(theta^2 * nk) / N - mu^2
  c(mu = mu, sigma = sqrt(max(sig2, 1e-4)))
}

# Fisher scoring for one freed item, with the smooth-L0 SBIC penalty
# when pen_logN > 0
em_mstep_item_g <- function(theta, njk_t, rjk_t, tps, free_tps,
                            a0, nu0, g0, est_a,
                            inner_maxit = 50L, inner_tol = 1e-9,
                            pen_logN = 0, pen_eps = NULL,
                            h0 = NULL, free_h_tps = integer(0),
                            lower_a = 0.1, upper_a = 10) {
  a <- a0; nu <- nu0; g <- g0
  h <- if (is.null(h0)) rep(0, length(g0)) else h0
  nfg <- length(free_tps)
  nfh <- length(free_h_tps)
  pen_on <- !is.null(pen_eps) && pen_logN > 0 && nfg > 0
  np <- as.integer(est_a) + 1L + nfg + nfh
  pos_a <- if (est_a) 1L else 0L
  pos_nu <- pos_a + 1L
  pos_h0 <- pos_nu + nfg
  q_obj <- function(a, nu, g, h) {
    val <- 0
    for (t in tps) {
      p <- pmin(pmax(stats::plogis((a + h[t]) * theta - nu - g[t]), 1e-10),
                1 - 1e-10)
      val <- val + sum(rjk_t[[t]] * log(p) +
                       (njk_t[[t]] - rjk_t[[t]]) * log(1 - p))
    }
    if (pen_on) {
      gf <- g[free_tps]
      val <- val - pen_logN * sum(gf^2 / (gf^2 + pen_eps))
    }
    val
  }
  newton_from <- function(a, nu, g, h, use_pen) {
    for (step_i in seq_len(inner_maxit)) {
      s <- numeric(np)
      I <- matrix(0, np, np)
      for (t in tps) {
        p <- pmin(pmax(stats::plogis((a + h[t]) * theta - nu - g[t]), 1e-10),
                  1 - 1e-10)
        resid <- rjk_t[[t]] - njk_t[[t]] * p
        w <- njk_t[[t]] * p * (1 - p)
        if (est_a) {
          s[pos_a] <- s[pos_a] + sum(resid * theta)
          I[pos_a, pos_a] <- I[pos_a, pos_a] + sum(w * theta^2)
          I[pos_a, pos_nu] <- I[pos_a, pos_nu] - sum(w * theta)
        }
        s[pos_nu] <- s[pos_nu] - sum(resid)
        I[pos_nu, pos_nu] <- I[pos_nu, pos_nu] + sum(w)
        gi <- match(t, free_tps)
        if (!is.na(gi)) {
          pg <- pos_nu + gi
          s[pg] <- -sum(resid)
          if (est_a) I[pos_a, pg] <- -sum(w * theta)
          I[pos_nu, pg] <- sum(w)
          I[pg, pg] <- sum(w)
          if (use_pen) {
            gt <- g[t]
            den <- (gt^2 + pen_eps)
            s[pg] <- s[pg] - pen_logN * 2 * gt * pen_eps / den^2
            d2 <- pen_logN * 2 * pen_eps * (pen_eps - 3 * gt^2) / den^3
            I[pg, pg] <- I[pg, pg] + max(d2, 0)
          }
        }
        hi <- match(t, free_h_tps)
        if (!is.na(hi)) {
          ph <- pos_h0 + hi
          s[ph] <- sum(resid * theta)
          I[ph, ph] <- sum(w * theta^2)
          if (est_a) I[pos_a, ph] <- sum(w * theta^2)
          I[pos_nu, ph] <- -sum(w * theta)
          if (!is.na(gi)) I[pos_nu + gi, ph] <- -sum(w * theta)
        }
      }
      I[lower.tri(I)] <- t(I)[lower.tri(I)]
      dp <- tryCatch(solve(I, s), error = function(e) NULL)
      if (is.null(dp)) break
      q0 <- if (use_pen) q_obj(a, nu, g, h) else NA_real_
      stepf <- 1
      accepted <- FALSE
      a1 <- a; nu1 <- nu; g1 <- g; h1 <- h
      for (hh in 1:12) {
        a1 <- if (est_a) a + stepf * dp[pos_a] else a
        nu1 <- nu + stepf * dp[pos_nu]
        g1 <- g
        h1 <- h
        if (nfg) g1[free_tps] <- g[free_tps] + stepf * dp[pos_nu + seq_len(nfg)]
        if (nfh) h1[free_h_tps] <- h[free_h_tps] +
            stepf * dp[pos_h0 + seq_len(nfh)]
        a1 <- min(max(a1, lower_a), upper_a)
        nu1 <- min(max(nu1, -10), 10)
        if (nfg) g1[free_tps] <- pmin(pmax(g1[free_tps], -10), 10)
        if (nfh) h1[free_h_tps] <- pmin(pmax(h1[free_h_tps], -10), 10)
        if (!use_pen || q_obj(a1, nu1, g1, h1) >= q0 - 1e-10) {
          accepted <- TRUE
          break
        }
        stepf <- stepf / 2
      }
      if (!accepted) break               # no improving step: keep old values
      chg <- max(abs(c(a1 - a, nu1 - nu,
                       if (nfg) g1[free_tps] - g[free_tps] else 0,
                       if (nfh) h1[free_h_tps] - h[free_h_tps] else 0)))
      a <- a1; nu <- nu1; g <- g1; h <- h1
      if (chg < inner_tol) break
    }
    list(a = a, nu = nu, g = g, h = h)
  }
  if (!pen_on) return(newton_from(a, nu, g, h, use_pen = FALSE))
  # two starts per update so items can switch between the shrunk
  # and the kept basin
  cand_a <- newton_from(a, nu, g, h, use_pen = TRUE)
  un <- newton_from(a, nu, g, h, use_pen = FALSE)
  cand_b <- newton_from(un$a, un$nu, un$g, un$h, use_pen = TRUE)
  if (q_obj(cand_b$a, cand_b$nu, cand_b$g, cand_b$h) >
      q_obj(cand_a$a, cand_a$nu, cand_a$g, cand_a$h)) cand_b else cand_a
}

# multiple group EM, group 1 fixed N(0,1), item parameters tied,
# freed items get per-group offsets
em_fit_multi <- function(data_list, model = c("2PL", "1PL"),
                         free_items = NULL, control = list(),
                         penalty = NULL, start = NULL,
                         free_slope_items = NULL, pweights = NULL) {
  model <- match.arg(model)
  ctrl <- em_control(control)
  md <- em_multi_data(data_list, collapse = ctrl$collapse,
                      pweights = pweights)
  J <- md$J
  n_groups <- md$n_groups
  reg <- em_model_registry(model)
  if (!is.null(ctrl$lower_a)) {
    reg$lower[["a"]] <- ctrl$lower_a
  }
  pm <- em_parmap_concurrent(md, model, free_items, free_slope_items)
  theta <- em_theta_grid(ctrl$theta_nodes, ctrl$theta_range)
  free_idx <- which(md$pool %in% union(free_items, free_slope_items))
  inv_idx <- setdiff(seq_len(J), free_idx)

  # start values from pooled item means across groups
  psum <- numeric(J); pn <- numeric(J)
  for (t in seq_len(n_groups)) {
    wv <- md$groups[[t]]
    psum[wv$idx] <- psum[wv$idx] + colSums(wv$ed$Y1 * wv$ed$w)
    pn[wv$idx] <- pn[wv$idx] + colSums(wv$ed$Obs * wv$ed$w)
  }
  pbar <- pmin(pmax(psum / pmax(pn, 1), 0.05), 0.95)
  a <- rep(1, J)
  nu <- -stats::qlogis(pbar)
  g <- matrix(0, J, n_groups)
  h <- matrix(0, J, n_groups)
  mu <- rep(0, n_groups); sig <- rep(1, n_groups)
  if (!is.null(start)) {
    a <- start$a; nu <- start$nu; g <- start$g
    if (!is.null(start$h)) h <- start$h
    mu <- start$mu; sig <- start$sigma
  }
  pis <- lapply(seq_len(n_groups), function(t)
    dnorm_discrete(theta, mu[t], sig[t]))
  # total smooth-L0 penalty, the loglik trace tracks loglik minus it
  pen_total <- function(g) {
    if (is.null(penalty) || length(free_idx) == 0) return(0)
    gf <- g[free_idx, -1, drop = FALSE]
    penalty$logN * sum(gf^2 / (gf^2 + penalty$eps))
  }
  est_a <- rep(reg$est[["a"]], J)
  lower_a <- rep(reg$lower[["a"]], J); upper_a <- rep(reg$upper[["a"]], J)
  lower_nu <- rep(reg$lower[["nu"]], J); upper_nu <- rep(reg$upper[["nu"]], J)

  ll_trace <- numeric(0)
  converged <- FALSE
  iter <- 0L
  for (s in seq_len(ctrl$maxit)) {
    es <- em_estep_multi(md, a, nu, g, theta, pis, h)
    obj <- es$loglik - pen_total(g)
    if (length(ll_trace) > 0 &&
        obj < ll_trace[length(ll_trace)] - 1e-8 * (1 + abs(obj))) {
      warning("EM log-likelihood decreased at iteration ", s, ".",
              call. = FALSE)
    }
    ll_trace <- c(ll_trace, obj)
    old <- c(a, nu, g[free_idx, , drop = FALSE],
             h[free_idx, , drop = FALSE], mu, sig)

    for (t in seq_len(n_groups)[-1]) {
      upd <- em_update_normal(es$groups[[t]]$nk, theta, md$groups[[t]]$N)
      mu[t] <- upd[["mu"]]; sig[t] <- upd[["sigma"]]
      pis[[t]] <- dnorm_discrete(theta, mu[t], sig[t])
    }

    if (length(inv_idx) > 0) {
      npool <- Reduce(`+`, lapply(es$groups, `[[`, "njk"))
      rpool <- Reduce(`+`, lapply(es$groups, `[[`, "rjk"))
      ms <- em_mstep_items_cpp(theta,
                               npool[inv_idx, , drop = FALSE],
                               rpool[inv_idx, , drop = FALSE],
                               a[inv_idx], nu[inv_idx], est_a[inv_idx],
                               lower_a[inv_idx], upper_a[inv_idx],
                               lower_nu[inv_idx], upper_nu[inv_idx],
                               ctrl$inner_maxit, ctrl$inner_tol)
      a[inv_idx] <- ms$a; nu[inv_idx] <- ms$nu
    }
    for (jj in free_idx) {
      it <- md$pool[jj]
      tps <- which(vapply(md$groups, function(w) jj %in% w$idx, logical(1)))
      # anchor (a, nu) at the first administration, otherwise
      # (nu, g) would be jointly unidentified
      free_tps <- if (it %in% free_items) tps[-1L] else integer(0)
      free_h_tps <- if (it %in% free_slope_items) tps[-1L]
                    else integer(0)
      njk_t <- lapply(es$groups, function(w) w$njk[jj, ])
      rjk_t <- lapply(es$groups, function(w) w$rjk[jj, ])
      res <- em_mstep_item_g(theta, njk_t, rjk_t, tps, free_tps,
                             a[jj], nu[jj], g[jj, ], est_a[jj],
                             inner_maxit = ctrl$inner_maxit,
                             inner_tol = ctrl$inner_tol,
                             pen_logN = if (is.null(penalty)) 0
                                        else penalty$logN,
                             pen_eps = if (is.null(penalty)) NULL
                                       else penalty$eps,
                             h0 = h[jj, ], free_h_tps = free_h_tps,
                             lower_a = lower_a[jj], upper_a = upper_a[jj])
      a[jj] <- res$a; nu[jj] <- res$nu; g[jj, ] <- res$g; h[jj, ] <- res$h
    }

    delta <- max(abs(c(a, nu, g[free_idx, , drop = FALSE],
                       h[free_idx, , drop = FALSE], mu, sig) - old))
    iter <- s
    if (ctrl$verbose)
      cat(sprintf("EM iter %d: loglik %.6f, max change %.2e\n",
                  s, es$loglik, delta))
    if (delta < ctrl$conv) { converged <- TRUE; break }
  }
  es <- em_estep_multi(md, a, nu, g, theta, pis, h)

  g_names <- character(0); g_vals <- numeric(0)
  h_names <- character(0); h_vals <- numeric(0)
  for (jj in free_idx) {
    it <- md$pool[jj]
    tps <- which(vapply(md$groups, function(w) jj %in% w$idx, logical(1)))
    for (t in tps[-1L]) {
      if (it %in% free_items) {
        g_names <- c(g_names, paste0(it, "__G", t, "_g"))
        g_vals <- c(g_vals, g[jj, t])
      }
      if (it %in% free_slope_items) {
        h_names <- c(h_names, paste0(it, "__G", t, "_h"))
        h_vals <- c(h_vals, h[jj, t])
      }
    }
  }
  npar <- sum(!duplicated(pm$parindex[pm$est])) + 2L * (n_groups - 1L)
  N_total <- sum(vapply(md$groups, `[[`, numeric(1), "N"))
  # under the smooth-L0 penalty the BIC counts the offsets by their
  # smooth indicator g^2 / (g^2 + eps), not as full parameters
  npar_bic <- npar
  if (!is.null(penalty)) {
    offs <- c(g_vals, h_vals)
    npar_bic <- npar - length(offs) +
      sum(offs^2 / (offs^2 + penalty$eps))
  }
  dats <- lapply(data_list, function(d) {
    X <- as.matrix(d); storage.mode(X) <- "integer"; X
  })
  structure(
    list(item = md$pool, a = as.numeric(a), nu = as.numeric(nu),
         b = as.numeric(nu / a), g = g, h = h,
         g_estimates = stats::setNames(g_vals, g_names),
         h_estimates = stats::setNames(h_vals, h_names),
         counts = es$groups,
         trend = data.frame(group = seq_len(n_groups), mu = mu, sigma = sig),
         model = model, parmap = pm, theta = theta, pis = pis,
         loglik = es$loglik, loglik_trace = ll_trace,
         deviance = -2 * es$loglik, npar = npar,
         bic = -2 * es$loglik + log(N_total) * npar_bic,
         converged = converged, iter = iter, control = ctrl,
         free_items = free_items, dat = dats, pweights = pweights),
    class = "irtlink_em"
  )
}
