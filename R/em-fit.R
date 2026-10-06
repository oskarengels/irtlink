# R/em-fit.R
# internal marginal maximum likelihood EM calibration engine
# (Bock & Aitkin, 1981), update formulas as in Hanson (1998)

# equally spaced theta grid, default 101 nodes on [-6, 6]
em_theta_grid <- function(nodes = 101L, range = c(-6, 6)) {
  seq(range[1], range[2], length.out = nodes)
}

# discrete normal density on a theta grid
dnorm_discrete <- function(x, mean, sd) {
  p <- stats::dnorm(x = x, mean = mean, sd = sd)
  p / sum(p)
}

# item probabilities plogis(a_j * theta_k - nu_j), clamped away
# from 0/1, with `guess` the 3PL form
em_irf_matrix <- function(a, nu, theta, guess = NULL) {
  P <- stats::plogis(outer(a, theta) - nu)
  if (!is.null(guess)) P <- guess + (1 - guess) * P
  pmin(pmax(P, 1e-10), 1 - 1e-10)
}

# control defaults for the EM loop
em_control <- function(control = list()) {
  ctrl <- list(theta_nodes = 101L, theta_range = c(-6, 6), maxit = 1000L,
               conv = 1e-5, inner_maxit = 5L, inner_tol = 1e-8,
               collapse = TRUE, verbose = FALSE, prior_c = c(5, 17),
               lower_a = NULL)
  unknown <- setdiff(names(control), names(ctrl))
  if (length(unknown) > 0) {
    stop("Unknown control option(s): ", paste(unknown, collapse = ", "),
         call. = FALSE)
  }
  ctrl <- utils::modifyList(ctrl, control)
  if (!is.null(ctrl$lower_a) &&
      (!is.numeric(ctrl$lower_a) || length(ctrl$lower_a) != 1L ||
       !is.finite(ctrl$lower_a) || ctrl$lower_a <= 0 || ctrl$lower_a > 1)) {
    stop("control$lower_a must be a single value in (0, 1].", call. = FALSE)
  }
  ctrl
}

# collapse identical response patterns into weighted unique rows,
# person weights are summed within a pattern
em_collapse_patterns <- function(X, w = NULL) {
  if (is.null(w)) w <- rep(1, nrow(X))
  key <- apply(X, 1L, paste, collapse = ",")
  first <- !duplicated(key)
  if (sum(first) >= 0.8 * nrow(X)) {
    return(list(X = X, w = w))
  }
  sums <- tapply(w, key, sum)
  list(X = X[first, , drop = FALSE],
       w = as.numeric(sums[key[first]]))
}

# precomputed response masks for the E step
em_estep_data <- function(X, w) {
  Obs <- (!is.na(X)) * 1
  Y1 <- X
  Y1[is.na(Y1)] <- 0L
  storage.mode(Y1) <- "double"
  Y1 <- Obs * Y1
  list(Y1 = Y1, Y0 = Obs - Y1, Obs = Obs, w = w, N = nrow(X))
}

# E step (Hanson, 1998, Eqs. 12-14) in matrix form
em_estep <- function(ed, logP, logQ, log_pi, return_posterior = FALSE) {
  LL <- ed$Y1 %*% logP + ed$Y0 %*% logQ
  LL <- sweep(LL, 2L, log_pi, "+")
  m <- LL[cbind(seq_len(ed$N), max.col(LL, ties.method = "first"))]
  W <- exp(LL - m)
  s <- rowSums(W)
  loglik <- sum(ed$w * (m + log(s)))
  W <- (ed$w / s) * W                      # row i scaled by w_i / s_i
  njk <- crossprod(ed$Obs, W)
  rjk <- crossprod(ed$Y1, W)
  dimnames(njk) <- dimnames(rjk) <- NULL
  out <- list(njk = njk, rjk = rjk, nk = colSums(W), loglik = loglik)
  if (return_posterior) out$posterior <- W / ed$w
  out
}

# start values, a = 1 and nu from the logit of the item means
em_start_values <- function(X) {
  pbar <- colMeans(X == 1L, na.rm = TRUE)
  pbar <- pmin(pmax(pbar, 0.05), 0.95)
  list(a = rep(1, ncol(X)), nu = -stats::qlogis(pbar))
}

# EM loop for one group, est_dist = TRUE frees the latent
# distribution and `fixed` holds the bank items (FIPC)
em_fit <- function(X, model = c("2PL", "1PL", "3PL"), control = list(),
                   fixed = NULL, est_dist = FALSE, pweights = NULL) {
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
  reg <- em_model_registry(model)
  if (!is.null(ctrl$lower_a)) {
    reg$lower[["a"]] <- ctrl$lower_a
  }
  if (model == "3PL" && (!is.null(fixed) || est_dist)) {
    stop("fixed items and est_dist are not yet supported for the 3PL.",
         call. = FALSE)
  }
  pm <- em_parmap_separate(colnames(X), model)
  theta <- em_theta_grid(ctrl$theta_nodes, ctrl$theta_range)
  sv <- em_start_values(X)
  a <- sv$a
  nu <- sv$nu
  gu <- if (model == "3PL") rep(unname(reg$start[["c"]]), J) else rep(0, J)
  # Bayes-modal objective, loglik plus beta log-prior on guessing
  pr_tot <- function(gu) {
    if (model != "3PL") return(0)
    sum((ctrl$prior_c[1] - 1) * log(gu) +
        (ctrl$prior_c[2] - 1) * log(1 - gu))
  }
  fix_idx <- integer(0)
  if (!is.null(fixed)) {
    fix_idx <- which(colnames(X) %in% fixed$item)
    mfix <- match(colnames(X)[fix_idx], fixed$item)
    a[fix_idx] <- fixed$a[mfix]
    nu[fix_idx] <- fixed$a[mfix] * fixed$b[mfix]
    pm$est[pm$item %in% fixed$item] <- FALSE
    pm$parindex[] <- NA_integer_
    pm$parindex[pm$est] <- seq_len(sum(pm$est))
  }
  free_j <- setdiff(seq_len(J), fix_idx)
  est_a <- rep(reg$est[["a"]], J)
  lower_a <- rep(reg$lower[["a"]], J); upper_a <- rep(reg$upper[["a"]], J)
  lower_nu <- rep(reg$lower[["nu"]], J); upper_nu <- rep(reg$upper[["nu"]], J)
  cp <- if (ctrl$collapse) em_collapse_patterns(X, pweights)
        else list(X = X, w = if (is.null(pweights)) rep(1, nrow(X))
                             else pweights)
  ed <- em_estep_data(cp$X, cp$w)
  N_eff <- sum(cp$w)
  mu <- 0; sig <- 1
  pi_k <- dnorm_discrete(theta, mean = mu, sd = sig)
  log_pi <- log(pi_k)
  ll_trace <- numeric(0)
  converged <- FALSE
  iter <- 0L
  for (s in seq_len(ctrl$maxit)) {
    P <- em_irf_matrix(a, nu, theta,
                       guess = if (model == "3PL") gu else NULL)
    es <- em_estep(ed, log(P), log(1 - P), log_pi)
    obj <- es$loglik + pr_tot(gu)
    if (length(ll_trace) > 0 &&
        obj < ll_trace[length(ll_trace)] - 1e-8 * (1 + abs(obj))) {
      warning("EM log-likelihood decreased at iteration ", s, ".",
              call. = FALSE)
    }
    ll_trace <- c(ll_trace, obj)
    old <- c(a[free_j], nu[free_j], gu[free_j], mu, sig)
    if (est_dist) {
      upd <- em_update_normal(es$nk, theta, N_eff)
      mu <- upd[["mu"]]; sig <- upd[["sigma"]]
      pi_k <- dnorm_discrete(theta, mean = mu, sd = sig)
      log_pi <- log(pi_k)
    }
    if (length(free_j) > 0) {
      if (model == "3PL") {
        for (j in free_j) {
          r3 <- em_mstep_item_3pl(theta, es$njk[j, ], es$rjk[j, ],
                                  a[j], nu[j], gu[j],
                                  prior_c = ctrl$prior_c,
                                  lower_a = lower_a[j],
                                  upper_a = upper_a[j])
          a[j] <- r3$a; nu[j] <- r3$nu; gu[j] <- r3$c
        }
      } else {
        ms <- em_mstep_items_cpp(theta,
                                 es$njk[free_j, , drop = FALSE],
                                 es$rjk[free_j, , drop = FALSE],
                                 a[free_j], nu[free_j], est_a[free_j],
                                 lower_a[free_j], upper_a[free_j],
                                 lower_nu[free_j], upper_nu[free_j],
                                 ctrl$inner_maxit, ctrl$inner_tol)
        a[free_j] <- ms$a
        nu[free_j] <- ms$nu
      }
    }
    delta <- max(abs(c(a[free_j], nu[free_j], gu[free_j], mu, sig) - old))
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
  P <- em_irf_matrix(a, nu, theta,
                     guess = if (model == "3PL") gu else NULL)
  es <- em_estep(ed, log(P), log(1 - P), log_pi)
  pm$value[pm$parname == "a"] <- a
  pm$value[pm$parname == "nu"] <- nu
  if (model == "3PL") pm$value[pm$parname == "c"] <- gu
  structure(
    list(item = colnames(X), a = as.numeric(a), nu = as.numeric(nu),
         b = as.numeric(nu / a), c = as.numeric(gu),
         model = model, parmap = pm,
         theta = theta, pi = pi_k,
         dist = c(mu = mu, sigma = sig),
         loglik = es$loglik,
         loglik_trace = ll_trace, deviance = -2 * es$loglik,
         npar = sum(pm$est) + 2L * est_dist,
         njk = es$njk, rjk = es$rjk, dat = X, pweights = pweights,
         converged = converged, iter = iter, control = ctrl),
    class = "irtlink_em"
  )
}
