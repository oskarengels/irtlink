# R/link_joint.R

link_joint <- function(x, method = c("haberman", "haebara", "sl",
                                     "phl",
                                     "haebara_pw", "sl_pw"),
                       pow = 2, use_intercepts = TRUE, eps = NULL,
                       theta = seq(-6, 6, length.out = 101),
                       weights = "normal_1",
                       item_weights = c("inverse_admin", "uniform"),
                       restricted = FALSE, anchor = NULL, ref = NULL, ...) {
  method <- match.arg(method)
  item_weights <- match.arg(item_weights)
  ipars <- as_ipars(x)
  all_groups <- sort(unique(ipars$group))
  if (!is.null(anchor)) {
    ipars <- ipars[ipars$item %in% anchor, , drop = FALSE]
    if (nrow(ipars) == 0) {
      stop("No anchor items found in the item parameters.", call. = FALSE)
    }
    dropped <- setdiff(all_groups, sort(unique(ipars$group)))
    if (length(dropped) > 0) {
      warning("Anchor filtering removed all items from group(s) ",
              paste(dropped, collapse = ", "),
              "; these groups are absent from the trend.", call. = FALSE)
    }
  }
  groups <- sort(unique(ipars$group))
  n_groups <- length(groups)
  if (n_groups < 2) {
    stop("Joint linking requires at least two groups.", call. = FALSE)
  }

  if (restricted) {
    trend <- data.frame(group = groups,
                        mu = c(0, rep(NA_real_, n_groups - 1)),
                        sigma = c(1, rep(NA_real_, n_groups - 1)))
    converged <- c(TRUE, rep(NA, n_groups - 1))
    for (t in 2:n_groups) {
      sub <- ipars[ipars$group <= groups[t], , drop = FALSE]
      est <- joint_engine(sub, method, pow, use_intercepts, eps, theta,
                          weights, item_weights, ...)
      trend$mu[t] <- est$mu[t]
      trend$sigma[t] <- est$sigma[t]
      converged[t] <- est$converged
    }
    joint_fit <- NULL
  } else {
    est <- joint_engine(ipars, method, pow, use_intercepts, eps, theta,
                        weights, item_weights, ...)
    trend <- data.frame(group = groups, mu = est$mu, sigma = est$sigma)
    converged <- est$converged
    joint_fit <- est$joint_fit
  }

  trend <- rebase_trend(trend, ref)

  structure(
    list(
      method = method,
      approach = if (restricted) "joint_restricted" else "joint",
      trend = trend,
      steps = NULL,
      converged = converged,
      ipars = ipars,
      anchor = anchor,
      ref = ref,
      joint_fit = joint_fit,
      control = c(list(pow = pow, use_intercepts = use_intercepts,
                       eps = eps, theta = theta, weights = weights,
                       item_weights = item_weights, ref = ref),
                  list(...)),
      le = NULL,
      call = match.call()
    ),
    class = "irtlink_link"
  )
}

# dispatch one joint estimation over all groups, group 1 fixed at 0/1
joint_engine <- function(ipars, method, pow, use_intercepts, eps, theta,
                         weights = "normal_1", item_weights = "inverse_admin",
                         ...) {
  if (method %in% c("haberman", "haebara", "sl") &&
      !pow %in% c(0, 0.25, 0.5, 1, 2)) {
    stop("`pow` must be one of 0, 0.25, 0.5, 1, 2.", call. = FALSE)
  }
  if (any(ipars$a <= 0)) {
    stop("Non-positive discrimination(s) in item parameters; ",
         "remove degenerate items before linking.", call. = FALSE)
  }
  if (method == "phl") {
    return(joint_phl(ipars, item_weights = item_weights))
  }
  # the simultaneous SL criterion with free item parameters is not
  # identified, so "sl" runs the pairwise criterion
  if (method %in% c("sl", "haebara_pw", "sl_pw")) {
    return(est_joint_pairwise(
      ipars, method = if (method == "haebara_pw") "Hae" else "SL",
      theta = theta, weights = weights))
  }
  if (method == "haebara") {
    # the simultaneous criterion re-estimates free item parameters and
    # is derived without guessing, the pairwise form supports the 3PL
    if (!is.null(ipars$c) && any(ipars$c > 0)) {
      stop("method = \"haebara\" with variant = \"simultaneous\" does ",
           "not support guessing parameters; use variant = \"pairwise\".",
           call. = FALSE)
    }
    if (ipars_poly(ipars)) {
      stop("method = \"haebara\" with variant = \"simultaneous\" does ",
           "not support polytomous items; use approach = \"chain\" or ",
           "method = \"haberman\".", call. = FALSE)
    }
    if (is.null(eps)) eps <- link_eps(pow)
    return(est_joint_haebara(ipars, theta = theta, weights = weights,
                             pow = pow, eps = eps))
  }
  if (is.null(eps)) eps <- link_eps(pow)
  est_joint_haberman(ipars, pow = pow,
                     use_intercepts = use_intercepts, eps = eps)
}
