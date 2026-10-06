# R/le-sl.R
# closed-form approximate jackknife for chain Stocking-Lord linking
# (Robitzsch, 2025, Foundations 2, Eqs. 30-38)

# delete-one one-step estimates for one step, type = "symm" adds
# the second alignment direction
le_sl_step_targets <- function(a0, b0, a1, b1, mu, sigma, theta, wgt,
                               type = c("asymm", "symm")) {
  type <- match.arg(type)
  I <- length(a0)
  th1 <- mu + sigma * theta
  C1 <- vapply(seq_len(I), function(i)
    stats::plogis(a0[i] * (th1 - b0[i])), numeric(length(theta)))
  C2 <- vapply(seq_len(I), function(i)
    stats::plogis(a1[i] * (theta - b1[i])), numeric(length(theta)))
  P1d <- vapply(seq_len(I), function(i)
    a0[i] * C1[, i] * (1 - C1[, i]), numeric(length(theta)))
  Dmu <- rowMeans(P1d)
  Dsi <- Dmu * theta
  if (type == "symm") {
    psi <- (theta - mu) / sigma
    C1f <- vapply(seq_len(I), function(i)
      stats::plogis(a0[i] * (theta - b0[i])), numeric(length(theta)))
    C2p <- vapply(seq_len(I), function(i)
      stats::plogis(a1[i] * (psi - b1[i])), numeric(length(theta)))
    W2 <- vapply(seq_len(I), function(i)
      a1[i] * C2p[, i] * (1 - C2p[, i]), numeric(length(theta)))
    D2mu <- rowMeans(W2) / sigma
    D2si <- rowMeans(W2) * (theta - mu) / sigma^2
  }
  out <- matrix(NA_real_, I, 2)
  for (i in seq_len(I)) {
    Z <- C1[, i] - C2[, i]
    ci <- c(sum(wgt * Z * Dmu), sum(wgt * Z * Dsi))
    y <- rowSums(P1d[, -i, drop = FALSE])
    Bi <- rbind(c(sum(wgt * y * Dmu), sum(wgt * theta * y * Dmu)),
                c(sum(wgt * y * Dsi), sum(wgt * theta * y * Dsi)))
    if (type == "symm") {
      Z2 <- C1f[, i] - C2p[, i]
      ci <- ci + c(sum(wgt * Z2 * D2mu), sum(wgt * Z2 * D2si))
      y2mu <- rowSums(W2[, -i, drop = FALSE]) / sigma
      y2si <- rowSums(W2[, -i, drop = FALSE]) * (theta - mu) / sigma^2
      Bi <- Bi + rbind(c(sum(wgt * y2mu * D2mu), sum(wgt * y2si * D2mu)),
                       c(sum(wgt * y2mu * D2si), sum(wgt * y2si * D2si)))
    }
    shift <- tryCatch(
      solve(Bi, ci),
      error = function(e) {
        stop("AJK leave-one-out matrix for a Stocking-Lord step is ",
             "singular; linking error cannot be estimated. Original ",
             "error: ", conditionMessage(e), call. = FALSE)
      }
    )
    out[i, ] <- c(mu, sigma) + shift
  }
  out
}

# per step the closed-form delete-one estimates, then the chain
# accumulation through le_chain_target()
le_ajk_sl_chain <- function(link, method = "ajk") {
  st <- le_hae_settings(link)
  if (!identical(st$pow, 2)) {
    stop("AJK linking error for chain Stocking-Lord links is currently ",
         "available for pow = 2 only (Robitzsch, 2025, Foundations 2).",
         call. = FALSE)
  }
  from_metric <- identical(st$theta_metric, "from") &&
    identical(st$type, "asymm")
  comp <- le_ipar_components(link, unit_mode = "chain")
  ref_idx <- le_ref_idx(link, comp$groups)
  P <- comp$Tn - 1L
  delta_hat <- as.vector(rbind(link$steps$sigma_step, link$steps$mu_step))
  n_units <- length(comp$units)
  deltas <- matrix(delta_hat, nrow = n_units, ncol = 2L * P, byrow = TRUE)
  for (step in seq_len(P)) {
    sp <- 2L * step - 1L
    mp <- sp + 1L
    com <- intersect(comp$items_t[[step]], comp$items_t[[step + 1L]])
    com <- intersect(com, comp$units)
    if (!length(com)) next
    ix <- match(com, comp$units)
    mu_s <- link$steps$mu_step[step]
    sig_s <- link$steps$sigma_step[step]
    if (from_metric) {
      # expand in the swapped-fit coordinates, then invert back
      tg <- le_sl_step_targets(
        a0 = comp$aL[[step + 1L]][com], b0 = comp$bL[[step + 1L]][com],
        a1 = comp$aL[[step]][com], b1 = comp$bL[[step]][com],
        mu = -mu_s / sig_s, sigma = 1 / sig_s,
        theta = st$theta, wgt = st$wgt
      )
      deltas[ix, sp] <- 1 / tg[, 2L]
      deltas[ix, mp] <- -tg[, 1L] / tg[, 2L]
    } else {
      tg <- le_sl_step_targets(
        a0 = comp$aL[[step]][com], b0 = comp$bL[[step]][com],
        a1 = comp$aL[[step + 1L]][com], b1 = comp$bL[[step + 1L]][com],
        mu = mu_s, sigma = sig_s,
        theta = st$theta, wgt = st$wgt, type = st$type
      )
      deltas[ix, mp] <- tg[, 1L]
      deltas[ix, sp] <- tg[, 2L]
    }
  }
  full_target <- le_rebase_target_vec(
    le_chain_target(delta_hat, type = "mm", Tn = comp$Tn), ref_idx)
  targets <- t(apply(deltas, 1L, function(d) {
    le_rebase_target_vec(le_chain_target(d, type = "mm", Tn = comp$Tn),
                         ref_idx)
  }))
  se <- sqrt(n_units / (n_units - 1) *
               colSums((targets - rep(full_target, each = n_units))^2))
  le_from_target_se(se, groups = comp$groups, method = method,
                    n_reps = n_units, ref_idx = ref_idx)
}
