# R/joint-haebara.R
# joint Haebara linking, one simultaneous item response function
# criterion over all groups

joint_hae_prepare <- function(ipars, theta, wgt, pow, eps) {
  groups <- sort(unique(ipars$group))
  pool <- sort(unique(ipars$item))
  NI <- length(pool)
  NS <- length(groups)
  aM <- matrix(NA_real_, NI, NS)
  bM <- matrix(NA_real_, NI, NS)
  for (s in seq_len(NS)) {
    d <- ipars[ipars$group == groups[s], , drop = FALSE]
    aM[match(d$item, pool), s] <- d$a
    bM[match(d$item, pool), s] <- d$b
  }
  list(aM = aM, bM = bM, theta = theta, q = wgt / sum(wgt),
       pow = pow, eps = eps, NI = NI, NS = NS, groups = groups,
       pool = pool, is_1pl = all(abs(ipars$a - 1) < 1e-12))
}

joint_hae_unpack <- function(par, prep) {
  NI <- prep$NI
  NS <- prep$NS
  if (prep$is_1pl) {
    list(a = rep(1, NI), b = par[seq_len(NI)],
         mu = c(0, par[NI + seq_len(NS - 1)]), sigma = rep(1, NS))
  } else {
    list(a = par[seq_len(NI)], b = par[NI + seq_len(NI)],
         mu = c(0, par[2 * NI + seq_len(NS - 1)]),
         sigma = c(1, par[2 * NI + NS - 1 + seq_len(NS - 1)]))
  }
}

joint_hae_value <- function(par, prep) {
  p <- joint_hae_unpack(par, prep)
  val <- 0
  for (s in seq_len(prep$NS)) {
    est <- which(!is.na(prep$aM[, s]))
    if (!length(est)) next
    ah <- prep$aM[est, s]
    bh <- prep$bM[est, s]
    P_obs <- stats::plogis(outer(ah, prep$theta) - ah * bh)
    P_exp <- stats::plogis(outer(p$a[est] * p$sigma[s], prep$theta) -
                             p$a[est] * (p$b[est] - p$mu[s]))
    e <- P_obs - P_exp
    r <- if (prep$pow == 2) e^2 else (e^2 + prep$eps)^(prep$pow / 2)
    val <- val + sum(r %*% prep$q)
  }
  val
}

# analytic gradient of joint_hae_value in the packing order of
# joint_hae_unpack
joint_hae_gradient <- function(par, prep) {
  p <- joint_hae_unpack(par, prep)
  g_a <- numeric(prep$NI)
  g_b <- numeric(prep$NI)
  g_mu <- numeric(prep$NS)
  g_sigma <- numeric(prep$NS)
  for (s in seq_len(prep$NS)) {
    est <- which(!is.na(prep$aM[, s]))
    if (!length(est)) next
    ah <- prep$aM[est, s]
    bh <- prep$bM[est, s]
    P_obs <- stats::plogis(outer(ah, prep$theta) - ah * bh)
    P_exp <- stats::plogis(outer(p$a[est] * p$sigma[s], prep$theta) -
                             p$a[est] * (p$b[est] - p$mu[s]))
    e <- P_obs - P_exp
    r1 <- if (prep$pow == 2) 2 * e
          else prep$pow * e * (e^2 + prep$eps)^(prep$pow / 2 - 1)
    core <- -r1 * P_exp * (1 - P_exp)
    Wq <- sweep(core, 2L, prep$q, "*")
    thM <- matrix(prep$theta, nrow = length(est),
                  ncol = length(prep$theta), byrow = TRUE)
    rs <- rowSums(Wq)
    g_a[est] <- g_a[est] +
      rowSums(Wq * (p$sigma[s] * thM - (p$b[est] - p$mu[s])))
    g_b[est] <- g_b[est] - rs * p$a[est]
    g_mu[s] <- sum(rs * p$a[est])
    g_sigma[s] <- sum(rowSums(Wq * thM) * p$a[est])
  }
  if (prep$is_1pl) c(g_b, g_mu[-1L])
  else c(g_a, g_b, g_mu[-1L], g_sigma[-1L])
}

# start values from the column-centered difficulties and means
joint_hae_start <- function(prep) {
  b_mean <- colMeans(prep$bM, na.rm = TRUE)
  mu0 <- -(b_mean - b_mean[1L])[-1L]
  b0 <- rowMeans(sweep(prep$bM, 2L, b_mean, "-"), na.rm = TRUE)
  if (prep$is_1pl) c(b0, mu0)
  else c(rep(1, prep$NI), b0, mu0, rep(1, prep$NS - 1L))
}

est_joint_haebara <- function(ipars, theta, weights = "normal_1",
                              pow = 2, eps = 1e-3, maxit = 1000L) {
  w <- make_link_weights(weights, theta)
  prep <- joint_hae_prepare(ipars, w$theta, w$wgt, pow, eps)
  opt <- stats::optim(joint_hae_start(prep), fn = joint_hae_value,
                      gr = joint_hae_gradient, prep = prep,
                      method = "BFGS",
                      control = list(maxit = maxit, reltol = 1e-12))
  p <- joint_hae_unpack(opt$par, prep)
  list(mu = p$mu, sigma = p$sigma, converged = opt$convergence == 0,
       joint_fit = list(par = opt$par, prep = prep, value = opt$value,
                        item = prep$pool, a = p$a, b = p$b))
}
