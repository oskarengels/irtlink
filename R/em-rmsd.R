# R/em-rmsd.R
# data-based RMSD and the bias-corrected estimator T4
# (Robitzsch, 2025, Foundations 36, Eqs. 11-15 and 30)

# posterior weight matrix of one group under the fitted model
em_rmsd_posterior <- function(X, P, pi_k, wts = NULL) {
  X <- as.matrix(X)
  ed <- em_estep_data(X, if (is.null(wts)) rep(1, nrow(X)) else wts)
  es <- em_estep(ed, log(P), log(1 - P), log(pi_k),
                 return_posterior = TRUE)
  list(W = es$posterior * ed$w, Y1 = ed$Y1, Obs = ed$Obs)
}

# T0 and T4 per item, the node variances and covariances come from
# the estimating equations of the observed response functions
em_rmsd_stats <- function(njk, rjk, P, w, W, Y1, Obs) {
  den <- pmax(njk, 1e-10)
  p_obs <- rjk / den
  e <- p_obs - P
  rmsd2 <- as.numeric((e^2) %*% w)
  t0 <- sqrt(rmsd2)
  W2 <- W^2
  S2x <- crossprod(Y1, W2)
  S2n <- crossprod(Obs, W2)
  V_diag <- (S2x * (1 - 2 * p_obs) + p_obs^2 * S2n) / den^2
  bias <- as.numeric(V_diag %*% w)
  t4 <- numeric(length(t0))
  for (i in seq_along(t0)) {
    if (t0[i] < 1e-8) next
    Z <- (W * (Y1[, i] * Obs[, i])) -
      sweep(W * Obs[, i], 2L, p_obs[i, ], "*")
    Vi <- crossprod(Z) / tcrossprod(den[i, ])
    eo <- e[i, ] * w
    Ci <- as.numeric(crossprod(eo, Vi %*% eo))
    t4[i] <- max(t0[i] - bias[i] / (2 * t0[i]) + Ci / (2 * t0[i]^3), 0)
  }
  list(rmsd = t0, rmsd_bc = t4)
}

# RMSD table from a fitted multiple group em object, NA where an
# item was not administered
em_rmsd_table <- function(fit, drop_items = character(0)) {
  n_groups <- length(fit$counts)
  theta <- fit$theta
  wcols <- paste0("g", seq_len(n_groups))
  tab <- data.frame(item = fit$item, stringsAsFactors = FALSE)
  for (t in seq_len(n_groups)) {
    cnt <- fit$counts[[t]]
    admin <- rowSums(cnt$njk) > 0
    X <- as.matrix(fit$dat[[t]])
    idx <- match(colnames(X), fit$item)
    P_g <- em_irf_matrix(fit$a[idx] + fit$h[idx, t],
                         fit$nu[idx] + fit$g[idx, t], theta)
    wts <- if (is.null(fit$pweights)) NULL else fit$pweights[[t]]
    po <- em_rmsd_posterior(X, P_g, fit$pis[[t]], wts)
    st <- em_rmsd_stats(cnt$njk[idx, , drop = FALSE],
                        cnt$rjk[idx, , drop = FALSE],
                        P_g, fit$pis[[t]], po$W, po$Y1, po$Obs)
    r <- rep(NA_real_, length(fit$item))
    rb <- rep(NA_real_, length(fit$item))
    r[idx] <- st$rmsd
    rb[idx] <- st$rmsd_bc
    r[!admin] <- NA_real_
    rb[!admin] <- NA_real_
    tab[[wcols[t]]] <- r
    tab[[paste0(wcols[t], "_bc")]] <- rb
  }
  n_present <- rowSums(!is.na(tab[wcols]))
  tab <- tab[n_present >= 2, , drop = FALSE]
  safe_max <- function(v) {
    v <- v[is.finite(v)]
    if (length(v)) max(v) else NA_real_
  }
  tab$max <- apply(tab[wcols], 1, safe_max)
  tab$max_bc <- apply(tab[paste0(wcols, "_bc")], 1, safe_max)
  if (length(drop_items)) tab <- tab[!tab$item %in% drop_items, , drop = FALSE]
  rownames(tab) <- NULL
  tab
}

# data-based RMSD from a full-invariance multiple group fit
em_rmsd_data <- function(calib, control = list()) {
  rpw <- responses_per_group(calib)
  fit <- em_fit_multi(rpw, model = calib$model, control = control)
  em_rmsd_table(fit)
}

# RMSD with split_items fully freed, only anchor candidates remain
em_rmsd_data_split <- function(data_list, split_items, model = "2PL",
                               control = list()) {
  fit <- em_fit_multi(data_list, model = model,
                      free_items = split_items,
                      free_slope_items = split_items, control = control)
  em_rmsd_table(fit, drop_items = split_items)
}
