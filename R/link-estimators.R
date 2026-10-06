# R/link-estimators.R

# one linking step from group `from` to group `to`, wide format with
# one row per common item

# robust eps default for Lp linking
link_eps <- function(pow) if (pow == 0) 1e-2 else 1e-3

# mean-mean, sigma from arithmetic means of the discriminations
est_link_mm <- function(pars, ...) {
  if (any(pars$a_from <= 0) || any(pars$a_to <= 0)) {
    stop("Non-positive discrimination(s) in linking pair; ",
         "remove degenerate items before linking.", call. = FALSE)
  }
  sigma <- mean(pars$a_to) / mean(pars$a_from)
  mu <- mean(pars$b_from) - sigma * mean(pars$b_to)
  list(mu = mu, sigma = sigma, converged = TRUE)
}

# mean-geometric mean, sigma from the geometric mean of a_to / a_from
est_link_mgm <- function(pars, ...) {
  if (any(pars$a_from <= 0) || any(pars$a_to <= 0)) {
    stop("Non-positive discrimination(s) in linking pair; ",
         "remove degenerate items before linking.", call. = FALSE)
  }
  sigma <- exp(mean(log(pars$a_to)) - mean(log(pars$a_from)))
  mu <- mean(pars$b_from) - sigma * mean(pars$b_to)
  list(mu = mu, sigma = sigma, converged = TRUE)
}

# Haberman step linking with Lq loss on the pair
est_link_haberman <- function(pars, pow = 2, use_intercepts = TRUE, eps = NULL, ...) {
  if (!pow %in% c(0, 0.25, 0.5, 1, 2)) {
    stop("`pow` must be one of 0, 0.25, 0.5, 1, 2.", call. = FALSE)
  }
  if (any(pars$a_from <= 0) || any(pars$a_to <= 0)) {
    stop("Non-positive discrimination(s) in linking pair; ",
         "remove degenerate items before linking.", call. = FALSE)
  }
  if (is.null(eps)) eps <- link_eps(pow)
  ipars <- data.frame(
    group = rep(1:2, each = nrow(pars)),
    item = rep(pars$item, 2),
    a = c(pars$a_from, pars$a_to),
    b = c(pars$b_from, pars$b_to)
  )
  res <- est_joint_haberman(ipars, pow = pow,
                            use_intercepts = use_intercepts, eps = eps)
  list(mu = res$mu[2L], sigma = res$sigma[2L], converged = res$converged)
}

# response function step linking (Haebara and Stocking-Lord) on the
# internal two-group fit
est_link_response <- function(pars, method, pow = 2, weights = "uniform",
                           theta = seq(-6, 6, length.out = 101), eps = NULL,
                           type = c("asymm", "symm"),
                           theta_metric = c("to", "from"), ...) {
  type <- match.arg(type)
  theta_metric <- match.arg(theta_metric)
  if (type == "symm" && theta_metric == "from") {
    stop("`theta_metric` applies to the asymmetric type only; the ",
         "symmetric criterion uses both metrics.", call. = FALSE)
  }
  if (any(pars$a_from <= 0) || any(pars$a_to <= 0)) {
    stop("Non-positive discrimination(s) in linking pair; ",
         "remove degenerate items before linking.", call. = FALSE)
  }
  if (!pow %in% c(0, 0.25, 0.5, 1, 2)) {
    stop("`pow` must be one of 0, 0.25, 0.5, 1, 2.", call. = FALSE)
  }
  w <- make_link_weights(weights, theta)
  swap <- type == "asymm" && theta_metric == "from"
  if (is.null(eps)) eps <- link_eps(pow)
  # theta_metric = "from" swaps the group roles and inverts the
  # fitted transformation
  tf_cols <- grep("^tau[0-9]+_from$", names(pars), value = TRUE)
  tt_cols <- grep("^tau[0-9]+_to$", names(pars), value = TRUE)
  tau_from <- as.matrix(pars[tf_cols])
  tau_to <- as.matrix(pars[tt_cols])
  if (length(tf_cols) > 0 && any(!is.na(tau_from))) {
    res <- if (swap) {
      fit_link_response_poly(a0 = pars$a_to, b0 = pars$b_to,
                             tau0 = tau_to,
                             a1 = pars$a_from, b1 = pars$b_from,
                             tau1 = tau_from,
                             method = method, type = type, pow = pow,
                             eps = eps, theta = w$theta, wgt = w$wgt)
    } else {
      fit_link_response_poly(a0 = pars$a_from, b0 = pars$b_from,
                             tau0 = tau_from,
                             a1 = pars$a_to, b1 = pars$b_to,
                             tau1 = tau_to,
                             method = method, type = type, pow = pow,
                             eps = eps, theta = w$theta, wgt = w$wgt)
    }
    mu <- res$mu
    sig <- res$sigma
    if (swap) {
      mu <- -mu / sig
      sig <- 1 / sig
    }
    return(list(mu = mu, sigma = sig, converged = res$converged,
                res = res))
  }
  c_from <- if (is.null(pars$c_from)) rep(0, nrow(pars)) else pars$c_from
  c_to <- if (is.null(pars$c_to)) rep(0, nrow(pars)) else pars$c_to
  res <- if (swap) {
    fit_link_response(a0 = pars$a_to, b0 = pars$b_to,
                   a1 = pars$a_from, b1 = pars$b_from,
                   c0 = c_to, c1 = c_from,
                   method = method, type = type, pow = pow, eps = eps,
                   theta = w$theta, wgt = w$wgt)
  } else {
    fit_link_response(a0 = pars$a_from, b0 = pars$b_from,
                   a1 = pars$a_to, b1 = pars$b_to,
                   c0 = c_from, c1 = c_to,
                   method = method, type = type, pow = pow, eps = eps,
                   theta = w$theta, wgt = w$wgt)
  }
  mu <- res$mu
  sig <- res$sigma
  if (swap) {
    mu <- -mu / sig
    sig <- 1 / sig
  }
  list(mu = mu, sigma = sig, converged = res$converged, res = res)
}

est_link_haebara <- function(pars, ...) est_link_response(pars, method = "Hae", ...)

est_link_sl <- function(pars, ...) est_link_response(pars, method = "SL", ...)
