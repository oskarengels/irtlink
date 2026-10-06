# R/em-engine.R
# wrappers connecting the internal EM engine to calibrate()

# fit one group under a fixed N(0,1) distribution (separate path)
fit_em_separate <- function(X, model = "2PL", group = 1L, verbose = FALSE,
                            control = list(), pweights = NULL) {
  control$verbose <- verbose
  fit <- if (model %in% c("GPCM", "PCM")) {
    em_fit_poly(X, model = model, control = control, pweights = pweights)
  } else {
    em_fit(X, model = model, control = control, pweights = pweights)
  }
  list(ipars = extract_em_ipars(fit, group), model = fit,
       converged = fit$converged)
}

# concurrent (multiple group) calibration, full invariance with
# free_items = NULL, partial invariance otherwise
fit_em_concurrent <- function(data_list, model = "2PL", free_items = NULL,
                              verbose = FALSE, control = list(),
                              pweights = NULL) {
  control$verbose <- verbose
  fit <- em_fit_multi(data_list, model = model, free_items = free_items,
                      control = control, pweights = pweights)
  list(trend = fit$trend,
       ipars = concurrent_ipars(fit, data_list),
       model = fit, converged = fit$converged,
       g = fit$g_estimates, bic = fit$bic)
}

# sequential FIPC, each later group fixes the bank items, frees its
# distribution, and adds newly seen items to the bank
fit_em_fipc <- function(data_list, model = "2PL", verbose = FALSE,
                        control = list(), pweights = NULL) {
  n_groups <- length(data_list)
  w1 <- fit_em_separate(data_list[[1]], model = model, group = 1L,
                        verbose = verbose, control = control,
                        pweights = if (is.null(pweights)) NULL
                                   else pweights[[1]])
  bank <- data.frame(item = w1$ipars$item, a = w1$ipars$a, b = w1$ipars$b,
                     stringsAsFactors = FALSE)
  trend <- data.frame(group = seq_len(n_groups), mu = NA_real_,
                      sigma = NA_real_)
  trend$mu[1] <- 0; trend$sigma[1] <- 1
  converged <- logical(n_groups); converged[1] <- w1$converged
  ipars_all <- w1$ipars
  control$verbose <- verbose
  for (t in seq_len(n_groups)[-1]) {
    Xt <- as.matrix(as.data.frame(data_list[[t]]))
    fixed_t <- bank[bank$item %in% colnames(Xt), , drop = FALSE]
    ft <- em_fit(Xt, model = model, control = control,
                 fixed = fixed_t, est_dist = TRUE,
                 pweights = if (is.null(pweights)) NULL
                            else pweights[[t]])
    trend$mu[t] <- unname(ft$dist[["mu"]])
    trend$sigma[t] <- unname(ft$dist[["sigma"]])
    converged[t] <- ft$converged
    ip_t <- extract_em_ipars(ft, group = t)
    ipars_all <- rbind(ipars_all, ip_t)
    new_items <- setdiff(colnames(Xt), bank$item)
    if (length(new_items) > 0) {
      bank <- rbind(bank,
                    ip_t[ip_t$item %in% new_items, c("item", "a", "b")])
    }
  }
  list(trend = trend, ipars = ipars_all, model = NULL,
       converged = all(converged))
}

# regularized calibration, eps grid warm-started, minimum BIC kept
fit_em_sbic <- function(data_list, model = "2PL",
                        eps = c(0.01, 0.001, 0.0001),
                        verbose = FALSE, control = list(),
                        pweights = NULL) {
  control$verbose <- verbose
  presence <- table(unlist(lapply(data_list,
                                  function(d) colnames(as.data.frame(d)))))
  common <- names(presence)[presence >= 2L]
  if (length(common) == 0L)
    stop("SBIC calibration: no common items found across groups.",
         call. = FALSE)
  N_total <- if (is.null(pweights)) {
    sum(vapply(data_list, function(d) nrow(as.data.frame(d)), integer(1)))
  } else {
    sum(unlist(pweights))
  }
  best <- NULL; best_bic <- Inf; chosen <- NA_real_
  # warm-start the first penalized fit from the full-invariance fit
  inv <- em_fit_multi(data_list, model = model, control = control,
                      pweights = pweights)
  warm <- list(a = inv$a, nu = inv$nu, g = inv$g, h = inv$h,
               mu = inv$trend$mu, sigma = inv$trend$sigma)
  for (e in eps) {
    # the penalty acts on the loglik, half the BIC factor on the
    # deviance scale
    fit <- em_fit_multi(data_list, model = model, free_items = common,
                        control = control,
                        penalty = list(logN = log(N_total) / 2, eps = e),
                        start = warm, pweights = pweights)
    warm <- list(a = fit$a, nu = fit$nu, g = fit$g, h = fit$h,
                 mu = fit$trend$mu, sigma = fit$trend$sigma)
    if (is.finite(fit$bic) && fit$bic < best_bic) {
      best <- fit; best_bic <- fit$bic; chosen <- e
    }
  }
  if (is.null(best)) {
    stop("SBIC calibration: no eps value produced a finite BIC; ",
         "check the data or the eps grid.", call. = FALSE)
  }
  list(trend = best$trend,
       ipars = concurrent_ipars(best, data_list),
       model = best, converged = best$converged,
       g = best$g_estimates, bic = best_bic, eps = chosen)
}

# item parameters from an irtlink_em fit, polytomous fits add the
# centered threshold columns tau1, tau2, ...
extract_em_ipars <- function(fit, group) {
  low <- if (!is.null(fit$control$lower_a)) fit$control$lower_a else 0.1
  stuck <- fit$a <= low * 1.01
  if (any(stuck)) {
    warning("Group ", group, ": ", sum(stuck),
            " item(s) sit at the lower discrimination bound (a <= ",
            low, "); b estimates may be unreliable.", call. = FALSE)
  }
  cc <- fit$c
  if (is.null(cc) || !length(cc)) cc <- rep(0, length(fit$a))
  out <- data.frame(group = group, item = fit$item, a = fit$a, b = fit$b,
                    c = cc, stringsAsFactors = FALSE)
  if (!is.null(fit$tau)) {
    tau <- fit$tau
    colnames(tau) <- paste0("tau", seq_len(ncol(tau)))
    out <- cbind(out, as.data.frame(tau))
  }
  out
}
