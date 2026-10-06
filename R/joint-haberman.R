# R/joint-haberman.R
# joint Haberman linking with the smoothed Lq loss (Haberman, 2009,
# Lq extension as in Robitzsch, 2020, Stats)

hab_lq_weights <- function(r, pow, eps) {
  if (pow == 2) return(rep(1, length(r)))
  if (pow == 0) return(1 / (r^2 + eps))
  (r^2 + eps)^(pow / 2 - 1)
}

# exact weighted fit of the two-way model with col_1 = 0, solved
# through the group-effect Schur complement
hab_solve_twoway <- function(y0, w, d) {
  NS <- ncol(y0)
  R_i <- rowSums(w)
  if (any(R_i <= 0)) {
    stop("Haberman linking: an item has no observed parameters.",
         call. = FALSE)
  }
  Y_i <- rowSums(w * y0)
  wd <- w * d
  M <- crossprod(wd, wd / R_i)
  D_g <- colSums(wd * d)
  rhs <- colSums(wd * y0) - colSums(wd * (Y_i / R_i))
  A <- diag(D_g, NS) - M
  col_eff <- rep(0, NS)
  sol <- tryCatch(
    solve(A[-1L, -1L, drop = FALSE], rhs[-1L]),
    error = function(e) {
      stop("Haberman linking design is not connected enough to ",
           "identify all group effects. Original error: ",
           conditionMessage(e), call. = FALSE)
    }
  )
  col_eff[-1L] <- sol
  row_eff <- (Y_i - as.vector(wd %*% col_eff)) / R_i
  list(row = row_eff, col = col_eff)
}

# reweighted least squares around the exact two-way solver
hab_fit_lq <- function(y, d, pow, eps, maxit = 500L, tol = 1e-12) {
  obs <- !is.na(y)
  y0 <- ifelse(obs, y, 0)
  w <- matrix(0, nrow(y), ncol(y))
  w[obs] <- 1
  fit <- hab_solve_twoway(y0, w, d)
  if (pow == 2) {
    return(list(row = fit$row, col = fit$col, converged = TRUE))
  }
  converged <- FALSE
  for (iter in seq_len(maxit)) {
    r <- y0 - fit$row - d * outer(rep(1, nrow(y)), fit$col)
    w <- matrix(0, nrow(y), ncol(y))
    w[obs] <- hab_lq_weights(r[obs], pow, eps)
    fit_new <- hab_solve_twoway(y0, w, d)
    delta <- max(abs(fit_new$row - fit$row), abs(fit_new$col - fit$col))
    fit <- fit_new
    if (delta < tol) {
      converged <- TRUE
      break
    }
  }
  list(row = fit$row, col = fit$col, converged = converged)
}

est_joint_haberman <- function(ipars, pow = 2, use_intercepts = TRUE,
                               eps = 1e-3, maxit = 500L) {
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
  is_1pl <- all(abs(ipars$a - 1) < 1e-12)
  if (is_1pl) {
    sigma <- rep(1, NS)
    a_i <- rep(1, NI)
    conv_a <- TRUE
  } else {
    fit_a <- hab_fit_lq(log(aM), d = 1, pow = pow, eps = eps,
                        maxit = maxit)
    sigma <- exp(fit_a$col)
    a_i <- exp(fit_a$row)
    conv_a <- fit_a$converged
  }
  if (use_intercepts) {
    fit_b <- hab_fit_lq(aM * bM, d = -a_i, pow = pow, eps = eps,
                        maxit = maxit)
    mu <- fit_b$col
    b_i <- fit_b$row / a_i
  } else {
    fit_b <- hab_fit_lq(-sweep(bM, 2L, sigma, "*"), d = 1, pow = pow,
                        eps = eps, maxit = maxit)
    mu <- fit_b$col
    b_i <- -fit_b$row
  }
  list(mu = mu, sigma = sigma, converged = conv_a && fit_b$converged,
       joint_fit = NULL, item = pool, a = a_i, b = b_i)
}
