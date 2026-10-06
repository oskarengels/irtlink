# R/em-person.R
# person parameter estimation, expected a posteriori (EAP) and
# weighted likelihood (WLE, Warm, 1989)

# EAP and posterior SD per person
em_eap_core <- function(X, a, nu, theta, pi) {
  ed <- em_estep_data(X, rep(1, nrow(X)))
  P <- em_irf_matrix(a, nu, theta)
  es <- em_estep(ed, log(P), log(1 - P), log(pi), return_posterior = TRUE)
  W <- es$posterior
  est <- as.numeric(W %*% theta)
  m2 <- as.numeric(W %*% theta^2)
  data.frame(est = est, se = sqrt(pmax(m2 - est^2, 0)))
}

# WLE via vectorized Fisher scoring with the Warm correction
# J / (2 I), SE = 1 / sqrt(I)
em_wle_core <- function(X, a, nu, maxit = 50L, tol = 1e-8) {
  N <- nrow(X); J <- ncol(X)
  Obs <- (!is.na(X)) * 1
  Y1 <- X
  Y1[is.na(Y1)] <- 0L
  storage.mode(Y1) <- "double"
  Y1 <- Obs * Y1
  A1 <- matrix(a, N, J, byrow = TRUE)
  A2 <- matrix(a^2, N, J, byrow = TRUE)
  A3 <- matrix(a^3, N, J, byrow = TRUE)
  theta <- numeric(N)
  info <- rep(1, N)
  for (it in seq_len(maxit)) {
    P <- plogis(sweep(outer(theta, a), 2, nu, "-"))
    PQ <- P * (1 - P)
    info <- pmax(rowSums(Obs * A2 * PQ), 1e-10)
    jterm <- rowSums(Obs * A3 * PQ * (1 - 2 * P))
    score <- rowSums(A1 * (Y1 - Obs * P)) + jterm / (2 * info)
    delta <- pmin(pmax(score / info, -1), 1)
    theta <- pmin(pmax(theta + delta, -10), 10)
    if (max(abs(delta)) < tol) break
  }
  data.frame(est = theta, se = 1 / sqrt(info))
}

person_scores <- function(calib, method = c("eap", "wle")) {
  method <- match.arg(method)
  if (!inherits(calib, "irtlink_calib") || !identical(calib$engine, "em"))
    stop("person_scores requires a calibration with engine = \"em\".",
         call. = FALSE)
  if (!calib$calibration %in% c("separate", "concurrent"))
    stop("person_scores currently supports calibration = \"separate\" ",
         "or \"concurrent\".", call. = FALSE)
  if (is.null(calib$models))
    stop("person_scores needs the fitted models; re-run calibrate() ",
         "with keep_models = TRUE.", call. = FALSE)
  if (identical(calib$model, "3PL"))
    stop("person_scores does not yet support the 3PL.", call. = FALSE)
  if (calib$model %in% c("GPCM", "PCM")) {
    if (method == "wle")
      stop("WLE for polytomous models is not yet available; use ",
           "method = \"eap\".", call. = FALSE)
    out <- lapply(seq_along(calib$models), function(t) {
      fit <- calib$models[[t]]
      sc <- em_eap_core_poly(fit$dat, fit$a, fit$d, fit$K, fit$theta,
                             fit$pi)
      data.frame(group = t, person = seq_len(nrow(fit$dat)), sc)
    })
    return(do.call(rbind, out))
  }
  if (calib$calibration == "separate") {
    out <- lapply(seq_along(calib$models), function(t) {
      fit <- calib$models[[t]]
      sc <- if (method == "eap") {
        em_eap_core(fit$dat, fit$a, fit$nu, fit$theta, fit$pi)
      } else {
        em_wle_core(fit$dat, fit$a, fit$nu)
      }
      data.frame(group = t, person = seq_len(nrow(fit$dat)), sc)
    })
    return(do.call(rbind, out))
  }
  if (calib$calibration == "concurrent") {
    fit <- calib$models[[1]]
    out <- lapply(seq_along(fit$dat), function(t) {
      X <- fit$dat[[t]]
      idx <- match(colnames(X), fit$item)
      a_t <- fit$a[idx] + fit$h[idx, t]
      nu_t <- fit$nu[idx] + fit$g[idx, t]
      sc <- if (method == "eap") {
        em_eap_core(X, a_t, nu_t, fit$theta, fit$pis[[t]])
      } else {
        em_wle_core(X, a_t, nu_t)
      }
      data.frame(group = t, person = seq_len(nrow(X)), sc)
    })
    return(do.call(rbind, out))
  }
  stop("person_scores currently supports calibration = \"separate\" ",
       "or \"concurrent\".", call. = FALSE)
}
