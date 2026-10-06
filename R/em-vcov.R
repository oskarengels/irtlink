# R/em-vcov.R
# observed information covariance of the item parameters based on the
# Fisher identity, delta method from (a, nu) to (a, b)

# analytic marginal score per item, stacked as (a1, nu1, a2, nu2, ...)
em_score_ipars <- function(fit, a, nu, gu = NULL, ed = NULL) {
  if (is.null(ed)) ed <- em_estep_data(fit$dat, em_fit_weights(fit))
  three_pl <- identical(fit$model, "3PL")
  P <- em_irf_matrix(a, nu, fit$theta, guess = if (three_pl) gu else NULL)
  es <- em_estep(ed, log(P), log(1 - P), log(fit$pi))
  if (!three_pl) {
    resid <- es$rjk - es$njk * P
    s_a <- as.numeric(resid %*% fit$theta)
    s_nu <- -rowSums(resid)
    return(if (fit$model == "2PL") as.vector(rbind(s_a, s_nu)) else s_nu)
  }
  # 3PL weight does not cancel against dP/deta, the beta prior on c
  # adds its gradient
  psi <- stats::plogis(outer(a, fit$theta) - nu)
  Wt <- es$rjk / P - (es$njk - es$rjk) / (1 - P)
  dPdeta <- (1 - gu) * psi * (1 - psi)
  th <- matrix(fit$theta, nrow = nrow(P), ncol = ncol(P), byrow = TRUE)
  s_a <- rowSums(Wt * dPdeta * th)
  s_nu <- -rowSums(Wt * dPdeta)
  pc <- fit$control$prior_c
  s_c <- rowSums(Wt * (1 - psi)) +
    (pc[1] - 1) / gu - (pc[2] - 1) / (1 - gu)
  as.vector(rbind(s_a, s_nu, s_c))
}

# person weights of a fit, ones when it was unweighted
em_fit_weights <- function(fit) {
  if (is.null(fit$pweights)) rep(1, nrow(fit$dat)) else fit$pweights
}

em_vcov_supported <- function(fit) {
  if (!inherits(fit, "irtlink_em"))
    stop("`fit` must be an irtlink_em object.", call. = FALSE)
  if (!identical(unname(fit$dist), c(0, 1)))
    stop("em_vcov_ipars currently supports separate fits with a fixed ",
         "N(0,1) distribution.", call. = FALSE)
  invisible(TRUE)
}

# observed information of the free coefficients, symmetrized
# central differences of the analytic score
em_obs_info_ipars <- function(fit, h = 1e-5) {
  em_vcov_supported(fit)
  ed <- em_estep_data(fit$dat, em_fit_weights(fit))
  J <- length(fit$item)
  three_pl <- fit$model == "3PL"
  two_pl <- fit$model == "2PL"
  score <- function(par) {
    if (three_pl) {
      em_score_ipars(fit, par[seq(1, 3 * J, 3)], par[seq(2, 3 * J, 3)],
                     gu = par[seq(3, 3 * J, 3)], ed = ed)
    } else if (two_pl) {
      em_score_ipars(fit, par[seq(1, 2 * J, 2)], par[seq(2, 2 * J, 2)],
                     ed = ed)
    } else {
      em_score_ipars(fit, fit$a, par, ed = ed)
    }
  }
  par0 <- if (three_pl) as.vector(rbind(fit$a, fit$nu, fit$c))
          else if (two_pl) as.vector(rbind(fit$a, fit$nu))
          else fit$nu
  p <- length(par0)
  Jm <- matrix(0, p, p)
  for (j in seq_len(p)) {
    step <- h * (1 + abs(par0[j]))
    pp <- par0; pp[j] <- pp[j] + step
    pm <- par0; pm[j] <- pm[j] - step
    Jm[, j] <- (score(pp) - score(pm)) / (2 * step)
  }
  -(Jm + t(Jm)) / 2
}

# delta method from (a, nu) to (a, b)
em_delta_ab <- function(fit) {
  J <- length(fit$item)
  if (identical(fit$model, "3PL")) {
    out_names <- as.vector(t(cbind(paste0(fit$item, "_a"),
                                   paste0(fit$item, "_b"),
                                   paste0(fit$item, "_c"))))
    D <- matrix(0, 3 * J, 3L * J, dimnames = list(out_names, NULL))
    for (jj in seq_len(J)) {
      ca <- 3 * jj - 2L; cv <- 3 * jj - 1L; cc <- 3 * jj
      D[ca, ca] <- 1
      D[cv, ca] <- -fit$b[jj] / fit$a[jj]
      D[cv, cv] <- 1 / fit$a[jj]
      D[cc, cc] <- 1
    }
    return(D)
  }
  two_pl <- fit$model == "2PL"
  p <- if (two_pl) 2L * J else J
  out_names <- as.vector(t(cbind(paste0(fit$item, "_a"),
                                 paste0(fit$item, "_b"))))
  D <- matrix(0, 2 * J, p, dimnames = list(out_names, NULL))
  for (jj in seq_len(J)) {
    ar <- 2 * jj - 1L; br <- 2 * jj
    if (two_pl) {
      ca <- 2 * jj - 1L; cv <- 2 * jj
      D[ar, ca] <- 1
      D[br, ca] <- -fit$b[jj] / fit$a[jj]
      D[br, cv] <- 1 / fit$a[jj]
    } else {
      D[br, jj] <- 1 / fit$a[jj]
    }
  }
  D
}

em_vcov_ipars <- function(fit, h = 1e-5) {
  if (!is.null(fit$model) && fit$model %in% c("GPCM", "PCM")) {
    stop("Item parameter covariances are not yet available for ",
         "polytomous models.", call. = FALSE)
  }
  I_obs <- em_obs_info_ipars(fit, h = h)
  V_coef <- pinv(I_obs)
  D <- em_delta_ab(fit)
  V <- D %*% V_coef %*% t(D)
  V <- (V + t(V)) / 2
  dimnames(V) <- list(rownames(D), rownames(D))
  V
}
