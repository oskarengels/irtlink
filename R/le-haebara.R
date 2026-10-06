# R/le-haebara.R
# item scores of the chain Haebara criterion for the approximate
# jackknife (Robitzsch, 2024, Stats, Eq. 14)

# quadrature and loss settings of a fitted Haebara chain link
le_hae_settings <- function(link) {
  # the item scores below are derived without guessing parameters
  # and for dichotomous response functions
  if (!is.null(link$ipars$c) && any(link$ipars$c > 0)) {
    stop("The approximate jackknife for Haebara and Stocking-Lord ",
         "linking does not support guessing parameters; use ",
         "estimator = \"jackknife\".", call. = FALSE)
  }
  if (!is.null(link$ipars) && ipars_poly(link$ipars)) {
    stop("The approximate jackknife for Haebara and Stocking-Lord ",
         "linking does not support polytomous items; use ",
         "estimator = \"jackknife\".", call. = FALSE)
  }
  ctrl <- link$control
  if (is.null(ctrl)) ctrl <- list()
  # exact names, $ would partial-match theta to theta_metric
  pow <- if (is.null(ctrl[["pow"]])) 2 else ctrl[["pow"]]
  weights <- if (is.null(ctrl[["weights"]])) "uniform"
             else ctrl[["weights"]]
  theta <- if (is.null(ctrl[["theta"]])) seq(-6, 6, length.out = 101)
           else ctrl[["theta"]]
  eps <- if (is.null(ctrl[["eps"]])) link_eps(pow) else ctrl[["eps"]]
  type <- if (is.null(ctrl[["type"]])) "asymm" else ctrl[["type"]]
  theta_metric <- if (is.null(ctrl[["theta_metric"]])) "to"
                  else ctrl[["theta_metric"]]
  w <- make_link_weights(weights, theta)
  list(pow = pow, eps = eps, theta = w$theta, wgt = w$wgt,
       type = type, theta_metric = theta_metric)
}

# step parameters in the coordinates the AJK expands in
le_hae_delta_variant <- function(link, st) {
  sig <- link$steps$sigma_step
  mu <- link$steps$mu_step
  if (identical(st$theta_metric, "from") && identical(st$type, "asymm")) {
    as.vector(rbind(log(1 / sig), -mu / sig))
  } else {
    as.vector(rbind(log(sig), mu))
  }
}

# invert each step back to the forward parameterization, then
# accumulate through le_chain_target()
le_hae_target_from <- function(dstar, Tn) {
  P <- Tn - 1L
  sstar <- dstar[2L * seq_len(P) - 1L]
  mstar <- dstar[2L * seq_len(P)]
  sigma_step <- exp(-sstar)
  mu_step <- -mstar * sigma_step
  le_chain_target(as.vector(rbind(sigma_step, mu_step)), type = "mm",
                  Tn = Tn)
}

# step parameters from link$steps keep the expansion at the
# fitted solution
le_hae_delta <- function(link) {
  as.vector(rbind(log(link$steps$sigma_step), link$steps$mu_step))
}

# item scores of the step estimating equations, one row per link
# unit (Robitzsch, 2024, Stats, Eq. 15)
le_hae_gmat <- function(delta, comp, st) {
  P <- comp$Tn - 1L
  M <- matrix(0, nrow = length(comp$units), ncol = 2L * P,
              dimnames = list(comp$units, NULL))
  swap <- identical(st$theta_metric, "from") && identical(st$type, "asymm")
  for (step in seq_len(P)) {
    sp <- 2L * step - 1L
    mp <- sp + 1L
    sigma <- exp(delta[sp])
    m <- delta[mp]
    com <- intersect(comp$items_t[[step]], comp$items_t[[step + 1L]])
    com <- intersect(com, comp$units)
    if (!length(com)) next
    if (swap) {
      a0 <- comp$aL[[step + 1L]][com]
      b0 <- comp$bL[[step + 1L]][com]
      a1 <- comp$aL[[step]][com]
      b1 <- comp$bL[[step]][com]
    } else {
      a0 <- comp$aL[[step]][com]
      b0 <- comp$bL[[step]][com]
      a1 <- comp$aL[[step + 1L]][com]
      b1 <- comp$bL[[step + 1L]][com]
    }
    ix <- match(com, comp$units)
    for (ii in seq_along(com)) {
      psi <- stats::plogis(a0[ii] * (sigma * st$theta + m - b0[ii]))
      P2 <- stats::plogis(a1[ii] * (st$theta - b1[ii]))
      dlt <- psi - P2
      rho1 <- st$pow * dlt * (dlt^2 + st$eps)^(st$pow / 2 - 1)
      dm <- a0[ii] * psi * (1 - psi)
      g_m <- sum(st$wgt * rho1 * dm)
      g_s <- sum(st$wgt * rho1 * dm * sigma * st$theta)
      if (identical(st$type, "symm")) {
        psi2 <- (st$theta - m) / sigma
        P1f <- stats::plogis(a0[ii] * (st$theta - b0[ii]))
        P2p <- stats::plogis(a1[ii] * (psi2 - b1[ii]))
        dlt2 <- P1f - P2p
        rho12 <- st$pow * dlt2 * (dlt2^2 + st$eps)^(st$pow / 2 - 1)
        w2 <- a1[ii] * P2p * (1 - P2p)
        g_m <- g_m + sum(st$wgt * rho12 * w2 / sigma)
        g_s <- g_s + sum(st$wgt * rho12 * w2 * (st$theta - m) / sigma)
      }
      M[ix[ii], mp] <- g_m
      M[ix[ii], sp] <- g_s
    }
  }
  M
}
