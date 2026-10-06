# R/linking_error.R

linking_error <- function(link,
                          method = c("jackknife", "jackknife_bc",
                                     "ajk",
                                     "ajk_bc", "sandwich_esw",
                                     "sandwich_osw", "sandwich_bosw",
                                     "sandwich_esw_bc", "sandwich_osw_bc",
                                     "sandwich_bosw_bc",
                                     "sandwich_osw_analytic"),
                          vcov = NULL, h = 1e-5,
                          jk_scale = c("finite_item", "classical"),
                          cluster = NULL, ...) {
  method <- match.arg(method)
  jk_scale <- match.arg(jk_scale)
  if (!inherits(link, "irtlink_link"))
    stop("`link` must be an irtlink_link object.", call. = FALSE)
  if (is.null(link$ipars))
    stop("Linking-error estimation needs item parameters in `link$ipars`.",
         call. = FALSE)
  if (!is.null(cluster)) {
    if (is.factor(cluster)) cluster <- stats::setNames(as.character(cluster),
                                                       names(cluster))
    if (!is.character(cluster) || is.null(names(cluster)) ||
        any(names(cluster) == "") || anyNA(cluster)) {
      stop("`cluster` must be a named character vector mapping item ",
           "names to testlet labels.", call. = FALSE)
    }
    if (!method %in% c("jackknife", "ajk")) {
      stop("Cluster-robust (testlet) linking errors are available with ",
           "method = \"jackknife\" or method = \"ajk\".", call. = FALSE)
    }
    ok_ajk <- (identical(link$approach, "chain") &&
                 link$method %in% c("mgm", "mm")) ||
      (identical(link$approach, "joint") && identical(link$method, "phl"))
    if (method == "ajk" && !ok_ajk) {
      stop("The approximate jackknife with `cluster` is available for ",
           "chain mgm/mm and joint PHL links; use method = ",
           "\"jackknife\" otherwise.", call. = FALSE)
    }
  }
  if (method != "jackknife" && jk_scale != "finite_item") {
    warning("`jk_scale` is only used by method = \"jackknife\" and is ",
            "ignored here.", call. = FALSE)
  }
  # the sandwich and bias-corrected estimators are derived for
  # dichotomous models
  if (ipars_poly(link$ipars) && !method %in% c("jackknife", "ajk")) {
    stop("Linking errors for polytomous links are currently available ",
         "with method = \"jackknife\" or method = \"ajk\".", call. = FALSE)
  }
  # a non-default reference group only changes the projection side,
  # see le_rebase_jac() and le_ref_target_fn()
  if (method == "jackknife") {
    link$le <- le_jackknife(link, method = method, jk_scale = jk_scale,
                            cluster = cluster, ...)
  } else if (method == "jackknife_bc") {
    dots <- list(...)
    if (length(dots) > 0) {
      warning("Additional arguments are ignored for bias-corrected ",
              "jackknife linking error.", call. = FALSE)
    }
    link$le <- le_jackknife_bc(link, kind = "jackknife", vcov = vcov,
                               h = h, method = method)
  } else if (method == "ajk") {
    dots <- list(...)
    if (length(dots) > 0) {
      warning("Additional arguments are ignored for AJK linking error.",
              call. = FALSE)
    }
    link$le <- le_ajk(link, method = method, cluster = cluster)
  } else if (method == "ajk_bc") {
    dots <- list(...)
    if (length(dots) > 0) {
      warning("Additional arguments are ignored for bias-corrected AJK ",
              "linking error.", call. = FALSE)
    }
    link$le <- le_jackknife_bc(link, kind = "ajk", vcov = vcov, h = h,
                               method = method)
  } else if (method == "sandwich_esw") {
    dots <- list(...)
    if (length(dots) > 0) {
      warning("Additional arguments are ignored for expected sandwich ",
              "linking error.", call. = FALSE)
    }
    link$le <- le_esw(link)
  } else if (method %in% c("sandwich_esw_bc", "sandwich_osw_bc",
                           "sandwich_bosw_bc")) {
    dots <- list(...)
    if (length(dots) > 0) {
      warning("Additional arguments are ignored for bias-corrected ",
              "sandwich linking error.", call. = FALSE)
    }
    link$le <- le_sandwich_bc(link, vcov = vcov, method = method)
  } else if (method == "sandwich_osw_analytic") {
    dots <- list(...)
    if (length(dots) > 0) {
      warning("Additional arguments are ignored for analytic observed ",
              "sandwich linking error.", call. = FALSE)
    }
    link$le <- le_osw_analytic(link, method = method)
  } else {
    dots <- list(...)
    if (length(dots) > 0) {
      warning("Additional arguments are ignored for sandwich linking error.",
              call. = FALSE)
    }
    link$le <- le_sandwich(link, method)
  }
  link
}

# delete-one link item jackknife for chain and joint links, with
# `cluster` the deleted units are whole testlets (e.g., Monseur &
# Berezner, 2007)
le_jackknife <- function(link, method = "jackknife",
                         jk_scale = "finite_item", cluster = NULL, ...) {
  items <- jackknife_link_items(link)
  units <- le_cluster_sets(items, cluster)
  if (length(units) < 2L)
    stop("Jackknife linking error needs at least two link units.",
         call. = FALSE)

  reps <- lapply(seq_along(units), function(u) {
    # a failing refit is treated like a non-finite replicate
    refit <- tryCatch(refit_link_without_item(link, units[[u]], ...),
                      error = function(e) NULL)
    if (is.null(refit)) {
      return(data.frame(item = names(units)[u],
                        group = link$trend$group,
                        mu = NA_real_, sigma = NA_real_,
                        row.names = NULL))
    }
    data.frame(item = names(units)[u], refit$trend, row.names = NULL)
  })
  reps <- do.call(rbind, reps)

  nonfinite <- !is.finite(reps$mu) | !is.finite(reps$sigma)
  if (any(nonfinite)) {
    warning(sum(nonfinite), " delete-one estimate(s) were not finite and ",
            "are discarded from the jackknife linking error (items: ",
            paste(unique(reps$item[nonfinite]), collapse = ", "), ").",
            call. = FALSE)
  }

  groups <- link$trend$group
  le <- lapply(groups, function(w) {
    sub <- reps[reps$group == w, , drop = FALSE]
    full <- link$trend[match(w, link$trend$group), , drop = FALSE]
    data.frame(
      group = w,
      le_mu = jackknife_sd(sub$mu, full = full$mu, scale = jk_scale),
      le_sigma = jackknife_sd(sub$sigma, full = full$sigma,
                              scale = jk_scale),
      method = method,
      n_reps = nrow(sub)
    )
  })
  le <- do.call(rbind, le)
  rownames(le) <- NULL
  le
}

jackknife_sd <- function(x, full = NULL, scale = "finite_item") {
  x <- x[is.finite(x)]
  n <- length(x)
  if (n < 2L) return(NA_real_)
  if (is.null(full) || !is.finite(full)) return(NA_real_)
  fac <- if (identical(scale, "classical")) (n - 1) / n else n / (n - 1)
  sqrt(fac * sum((x - full)^2))
}

# named list of item sets, one per testlet, unmapped items stay
# singletons
le_cluster_sets <- function(items, cluster = NULL) {
  if (is.null(cluster)) {
    return(stats::setNames(as.list(items), items))
  }
  lab <- unname(cluster[items])
  lab[is.na(lab)] <- items[is.na(lab)]
  split(items, factor(lab, levels = unique(lab)))
}

jackknife_link_items <- function(link) {
  ip <- as_ipars(link$ipars)
  if (!is.null(link$anchor)) {
    return(intersect(link$anchor, unique(ip$item)))
  }
  presence <- table(ip$item)
  names(presence)[presence >= 2L]
}

refit_link_without_item <- function(link, item, ...) {
  ip <- link$ipars[!link$ipars$item %in% item, , drop = FALSE]
  anchor <- if (is.null(link$anchor)) NULL else setdiff(link$anchor, item)
  ctrl <- link_control(link, list(...))
  if (identical(link$approach, "chain")) {
    return(do.call(link_chain, c(list(x = ip, method = link$method,
                                     anchor = anchor), ctrl)))
  }
  if (link$approach %in% c("joint", "joint_restricted")) {
    return(do.call(link_joint, c(list(x = ip, method = link$method,
                                     restricted = link$approach == "joint_restricted",
                                     anchor = anchor), ctrl)))
  }
  stop("Unsupported linking approach for jackknife linking error: ",
       link$approach, call. = FALSE)
}

link_control <- function(link, overrides) {
  ctrl <- link$control
  if (is.null(ctrl)) ctrl <- list()
  utils::modifyList(ctrl, overrides)
}

# approximate one-step jackknife from the full estimating equations
le_ajk <- function(link, method = "ajk", cluster = NULL) {
  if (identical(link$approach, "chain")) {
    if (!link$method %in% c("mgm", "mm", "haebara", "sl")) {
      stop("AJK linking error for chain links is currently available only ",
           "for method = \"mgm\", \"mm\", \"haebara\", or \"sl\".",
           call. = FALSE)
    }
    if (identical(link$method, "sl")) {
      # closed form (Robitzsch, 2025, Foundations 2)
      return(le_ajk_sl_chain(link, method = method))
    }
    # analytic derivatives for mgm/mm, numeric for haebara
    return(le_ajk_chain(link, method = method,
                        analytic = !identical(link$method, "haebara"),
                        cluster = cluster))
  }
  if (identical(link$approach, "joint")) {
    if (!identical(link$method, "phl")) {
      stop("AJK linking error for joint links is currently available only ",
           "for method = \"phl\".", call. = FALSE)
    }
    return(le_ajk_phl(link, method = method, analytic = TRUE,
                      cluster = cluster))
  }
  stop("AJK linking error is not currently available for approach = \"",
       link$approach, "\".", call. = FALSE)
}

# index sets of comp$units per testlet
le_cluster_idx <- function(units, cluster = NULL) {
  if (is.null(cluster)) return(NULL)
  sets <- le_cluster_sets(units, cluster)
  lapply(sets, function(it) match(it, units))
}

le_ajk_chain <- function(link, method = "ajk", analytic = FALSE,
                         cluster = NULL) {
  comp <- le_ipar_components(link, unit_mode = "chain")
  ref_idx <- le_ref_idx(link, comp$groups)
  if (identical(link$method, "haebara")) {
    if (isTRUE(analytic)) {
      stop("Analytic AJK derivatives are not available for chain Haebara ",
           "links.", call. = FALSE)
    }
    st <- le_hae_settings(link)
    from_metric <- identical(st$theta_metric, "from") &&
      identical(st$type, "asymm")
    target_fn <- if (from_metric) {
      function(d) le_hae_target_from(d, Tn = comp$Tn)
    } else {
      function(d) le_chain_target(d, type = "mgm", Tn = comp$Tn)
    }
    return(le_ajk_core(
      delta_hat = le_hae_delta_variant(link, st),
      gmat_fn = function(d) le_hae_gmat(d, comp, st),
      target_fn = le_ref_target_fn(target_fn, ref_idx),
      groups = comp$groups,
      units = comp$units,
      method = method,
      ref_idx = ref_idx
    ))
  }
  delta <- le_chain_delta(comp, type = link$method)
  if (isTRUE(analytic)) {
    return(le_ajk_core_analytic(
      delta_hat = delta,
      Gi = le_chain_gmat(delta, comp, type = link$method),
      A_full = le_chain_bread(delta, comp, type = link$method),
      unit_bread_fn = function(unit) {
        le_chain_bread_unit(delta, comp, type = link$method, unit = unit)
      },
      target_fn = le_ref_target_fn(
        function(d) le_chain_target(d, type = link$method, Tn = comp$Tn),
        ref_idx),
      groups = comp$groups,
      units = comp$units,
      method = method,
      ref_idx = ref_idx,
      unit_sets = le_cluster_idx(comp$units, cluster)
    ))
  }
  le_ajk_core(
    delta_hat = delta,
    gmat_fn = function(d) le_chain_gmat(d, comp, type = link$method),
    target_fn = le_ref_target_fn(
      function(d) le_chain_target(d, type = link$method, Tn = comp$Tn),
      ref_idx),
    groups = comp$groups,
    units = comp$units,
    method = method,
    ref_idx = ref_idx,
    unit_sets = le_cluster_idx(comp$units, cluster)
  )
}

le_ajk_phl <- function(link, method = "ajk", analytic = FALSE,
                       cluster = NULL) {
  comp <- le_ipar_components(link, unit_mode = "joint")
  ref_idx <- le_ref_idx(link, comp$groups)
  omega <- le_phl_weights(comp, variant = le_phl_variant(link))
  delta <- le_phl_delta(comp, omega)
  if (isTRUE(analytic)) {
    return(le_ajk_core_analytic(
      delta_hat = delta,
      Gi = le_phl_gmat(delta, comp, omega),
      A_full = le_phl_bread(delta, comp, omega),
      unit_bread_fn = function(unit) {
        le_phl_bread_unit(delta, comp, omega, unit = unit)
      },
      target_fn = le_ref_target_fn(
        function(d) le_phl_target(d, Tn = comp$Tn), ref_idx),
      groups = comp$groups,
      units = comp$units,
      method = method,
      ref_idx = ref_idx,
      unit_sets = le_cluster_idx(comp$units, cluster)
    ))
  }
  le_ajk_core(
    delta_hat = delta,
    gmat_fn = function(d) le_phl_gmat(d, comp, omega),
    target_fn = le_ref_target_fn(
      function(d) le_phl_target(d, Tn = comp$Tn), ref_idx),
    groups = comp$groups,
    units = comp$units,
    method = method,
    ref_idx = ref_idx,
    unit_sets = le_cluster_idx(comp$units, cluster)
  )
}

le_ajk_core_analytic <- function(delta_hat, Gi, A_full, unit_bread_fn,
                                 target_fn, groups, units, method = "ajk",
                                 ref_idx = NULL, unit_sets = NULL) {
  if (is.null(unit_sets)) unit_sets <- as.list(seq_along(units))
  targets <- le_ajk_targets_analytic(delta_hat, Gi, A_full, unit_bread_fn,
                                     target_fn, units, unit_sets)
  full_target <- target_fn(delta_hat)
  n_units <- length(unit_sets)
  se <- sqrt(n_units / (n_units - 1) *
               colSums((targets - rep(full_target, each = n_units))^2))
  le_from_target_se(se, groups = groups, method = method, n_reps = n_units,
                    ref_idx = ref_idx)
}

le_ajk_targets_analytic <- function(delta_hat, Gi, A_full, unit_bread_fn,
                                    target_fn, units, unit_sets) {
  full_target <- target_fn(delta_hat)
  n_units <- length(unit_sets)
  targets <- matrix(NA_real_, nrow = n_units, ncol = length(full_target))

  for (i in seq_len(n_units)) {
    idx <- unit_sets[[i]]
    H_unit <- Reduce(`+`, lapply(units[idx], unit_bread_fn))
    A_loo <- A_full - H_unit
    g_unit <- colSums(Gi[idx, , drop = FALSE])
    shift <- tryCatch(
      solve(A_loo, as.numeric(g_unit)),
      error = function(e) {
        stop("AJK leave-one-out bread matrix is singular; linking error ",
             "cannot be estimated. Original error: ", conditionMessage(e),
             call. = FALSE)
      }
    )
    targets[i, ] <- target_fn(delta_hat + as.numeric(shift))
  }
  targets
}

le_ajk_core <- function(delta_hat, gmat_fn, target_fn, groups, units,
                        method = "ajk", ref_idx = NULL, unit_sets = NULL) {
  if (is.null(unit_sets)) unit_sets <- as.list(seq_along(units))
  targets <- le_ajk_targets(delta_hat, gmat_fn, target_fn, unit_sets)
  full_target <- target_fn(delta_hat)
  n_units <- length(unit_sets)
  se <- sqrt(n_units / (n_units - 1) *
               colSums((targets - rep(full_target, each = n_units))^2))
  le_from_target_se(se, groups = groups, method = method, n_reps = n_units,
                    ref_idx = ref_idx)
}

le_ajk_targets <- function(delta_hat, gmat_fn, target_fn, unit_sets) {
  Gsum <- function(d) colSums(gmat_fn(d))
  A_full <- le_numjac(Gsum, delta_hat)
  Gi <- gmat_fn(delta_hat)
  n_units <- length(unit_sets)
  targets <- matrix(NA_real_, nrow = n_units, ncol = length(target_fn(delta_hat)))

  for (i in seq_len(n_units)) {
    idx <- unit_sets[[i]]
    Hi <- le_numjac(function(d) colSums(gmat_fn(d)[idx, , drop = FALSE]),
                    delta_hat)
    A_loo <- A_full - Hi
    g_unit <- colSums(Gi[idx, , drop = FALSE])
    shift <- tryCatch(
      solve(A_loo, as.numeric(g_unit)),
      error = function(e) {
        stop("AJK leave-one-out bread matrix is singular; linking error ",
             "cannot be estimated. Original error: ", conditionMessage(e),
             call. = FALSE)
      }
    )
    # deleting a unit leaves -sum of its g_i, the one-step shift is
    # A^-1 times that sum
    targets[i, ] <- target_fn(delta_hat + as.numeric(shift))
  }
  targets
}

# finite-N bias correction for the delete-one and approximate
# jackknife (Robitzsch, 2025, Foundations 2, Eqs. 38-44)
le_jackknife_bc <- function(link, kind = c("jackknife", "ajk"), vcov,
                            h = 1e-5, method = "jackknife_bc",
                            derivatives = c("analytic", "numeric")) {
  kind <- match.arg(kind)
  derivatives <- match.arg(derivatives)
  if (is.null(vcov)) {
    stop("Bias-corrected jackknife linking error needs `vcov`, the ",
         "item-parameter covariance matrix/list.", call. = FALSE)
  }
  ctx <- le_bc_context(link)
  gamma <- le_gamma_vector(ctx$comp)
  Vgamma <- le_vcov_matrix(vcov, ctx$comp)

  full_fun <- function(g) le_target_from_gamma(g, ctx)
  target0 <- full_fun(gamma)
  loo_targets <- le_loo_targets_from_gamma(gamma, ctx, kind = kind)
  n_units <- nrow(loo_targets)
  diffs <- loo_targets - rep(target0, each = n_units)
  V_LE <- n_units / (n_units - 1) * crossprod(diffs)

  V_bias <- matrix(0, nrow = length(target0), ncol = length(target0))
  if (derivatives == "analytic") {
    full_ip_key <- paste(ctx$comp$ip$group, ctx$comp$ip$item)
    n_gamma <- 2L * nrow(ctx$comp$ip)
    U_full <- le_bc_U_analytic(ctx, ctx$comp, full_ip_key, n_gamma)
    V_SE <- U_full %*% Vgamma %*% t(U_full)
    for (i in seq_along(ctx$comp$units)) {
      unit <- ctx$comp$units[i]
      comp_i <- le_comp_drop_unit(ctx$comp, unit)
      U_i <- le_bc_U_analytic(ctx, comp_i, full_ip_key, n_gamma)
      U_diff <- U_i - U_full
      V_bias <- V_bias + U_diff %*% Vgamma %*% t(U_diff)
    }
  } else {
    U_full <- le_numjac(full_fun, gamma, h = h)
    V_SE <- U_full %*% Vgamma %*% t(U_full)
    for (i in seq_along(ctx$comp$units)) {
      unit <- ctx$comp$units[i]
      U_i <- le_numjac(function(g) {
        le_loo_target_from_gamma(g, ctx, kind = "jackknife", unit = unit)
      }, gamma, h = h)
      U_diff <- U_i - U_full
      V_bias <- V_bias + U_diff %*% Vgamma %*% t(U_diff)
    }
  }
  V_bias <- n_units / (n_units - 1) * V_bias
  V_LEbc <- V_LE - V_bias

  le_error_table(V_SE = V_SE, V_LE = V_LE, V_LEbc = V_LEbc,
                 groups = ctx$comp$groups, method = method,
                 n_reps = n_units, ref_idx = ctx$ref_idx)
}

# analytic U matrices for the bias correction,
# d target / d gamma = -J_T A^{-1} H_gamma
# (Robitzsch, 2025, Foundations 2, Eqs. 40-41)
le_bc_U_analytic <- function(ctx, comp, full_ip_key, n_gamma) {
  spec <- le_link_spec(comp, ctx)
  delta <- spec$delta
  if (ctx$approach == "chain") {
    A <- le_chain_bread(delta, comp, type = ctx$method)
    H <- le_chain_hgamma(delta, comp, type = ctx$method,
                         full_ip_key = full_ip_key, n_gamma = n_gamma)
    Jt <- le_rebase_jac(le_chain_target(delta, type = ctx$method,
                                        Tn = comp$Tn), ctx$ref_idx) %*%
      le_chain_target_jac(delta, type = ctx$method, Tn = comp$Tn)
  } else {
    omega <- le_phl_weights(comp, variant = ctx$phl_weights)
    A <- le_phl_bread(delta, comp, omega)
    H <- le_phl_hgamma(delta, comp, omega,
                       full_ip_key = full_ip_key, n_gamma = n_gamma)
    Jt <- le_rebase_jac(le_phl_target(delta, Tn = comp$Tn),
                        ctx$ref_idx) %*%
      le_phl_target_jac(delta, Tn = comp$Tn)
  }
  shift <- tryCatch(
    solve(A, H),
    error = function(e) {
      stop("Bias-correction bread matrix is singular; linking error ",
           "cannot be estimated. Original error: ", conditionMessage(e),
           call. = FALSE)
    }
  )
  -Jt %*% shift
}

# d G / d gamma for chain MGM/MM, only the item's own parameters
# at the two step groups carry derivatives
le_chain_hgamma <- function(delta, comp, type = c("mgm", "mm"),
                            full_ip_key, n_gamma) {
  type <- match.arg(type)
  P <- comp$Tn - 1L
  H <- matrix(0, nrow = 2L * P, ncol = n_gamma)
  for (step in seq_len(P)) {
    sp <- 2L * step - 1L
    mp <- sp + 1L
    com <- intersect(comp$items_t[[step]], comp$items_t[[step + 1L]])
    com <- intersect(com, comp$units)
    if (!length(com)) next
    r0 <- match(paste(comp$groups[step], com), full_ip_key)
    r1 <- match(paste(comp$groups[step + 1L], com), full_ip_key)
    a0 <- comp$aL[[step]][com]
    a1 <- comp$aL[[step + 1L]][com]
    if (type == "mgm") {
      H[sp, 2L * r0 - 1L] <- H[sp, 2L * r0 - 1L] - 1 / a0
      H[sp, 2L * r1 - 1L] <- H[sp, 2L * r1 - 1L] + 1 / a1
      H[mp, 2L * r0] <- H[mp, 2L * r0] - 1
      H[mp, 2L * r1] <- H[mp, 2L * r1] + exp(delta[sp])
    } else {
      H[sp, 2L * r0 - 1L] <- H[sp, 2L * r0 - 1L] + delta[sp]
      H[sp, 2L * r1 - 1L] <- H[sp, 2L * r1 - 1L] - 1
      H[mp, 2L * r0] <- H[mp, 2L * r0] - 1
      H[mp, 2L * r1] <- H[mp, 2L * r1] + delta[sp]
    }
  }
  H
}

# d G / d gamma for joint PHL, mirrors le_phl_gmat
le_phl_hgamma <- function(delta, comp, omega, full_ip_key, n_gamma) {
  P <- comp$Tn - 1L
  s <- c(0, delta[seq_len(P)])
  sigma <- exp(s)
  H <- matrix(0, nrow = 2L * P, ncol = n_gamma)
  for (ii in seq_along(comp$units)) {
    it <- comp$units[ii]
    tps <- which(vapply(comp$items_t, function(v) it %in% v, logical(1)))
    if (length(tps) < 2L) next
    w <- omega[[it]]
    rr <- match(paste(comp$groups[tps], it), full_ip_key)
    for (k in tps[tps >= 2L]) {
      oth <- setdiff(tps, k)
      r_s <- k - 1L
      r_m <- P + k - 1L
      kpos <- match(k, tps)
      H[r_s, 2L * rr[kpos] - 1L] <- H[r_s, 2L * rr[kpos] - 1L] -
        w * length(oth) / comp$aL[[k]][[it]]
      H[r_m, 2L * rr[kpos]] <- H[r_m, 2L * rr[kpos]] +
        w * length(oth) * sigma[k]
      for (hh in oth) {
        hpos <- match(hh, tps)
        H[r_s, 2L * rr[hpos] - 1L] <- H[r_s, 2L * rr[hpos] - 1L] +
          w / comp$aL[[hh]][[it]]
        H[r_m, 2L * rr[hpos]] <- H[r_m, 2L * rr[hpos]] - w * sigma[hh]
      }
    }
  }
  H
}

le_bc_context <- function(link,
                          what = "Bias-corrected jackknife linking error") {
  if (identical(link$approach, "chain")) {
    if (!link$method %in% c("mgm", "mm")) {
      stop(what, " for chain links is ",
           "currently available only for method = \"mgm\" or \"mm\".",
           call. = FALSE)
    }
    comp <- le_ipar_components(link, unit_mode = "chain")
    return(list(comp = comp, approach = "chain", method = link$method,
                ref_idx = le_ref_idx(link, comp$groups)))
  }
  if (identical(link$approach, "joint")) {
    if (!identical(link$method, "phl")) {
      stop(what, " for joint links is ",
           "currently available only for method = \"phl\".",
           call. = FALSE)
    }
    comp <- le_ipar_components(link, unit_mode = "joint")
    return(list(comp = comp, approach = "joint", method = link$method,
                phl_weights = le_phl_variant(link),
                ref_idx = le_ref_idx(link, comp$groups)))
  }
  stop(what, " is not currently available ",
       "for approach = \"", link$approach, "\".", call. = FALSE)
}

le_loo_targets_from_gamma <- function(gamma, ctx,
                                      kind = c("jackknife", "ajk")) {
  kind <- match.arg(kind)
  if (kind == "ajk") {
    return(le_ajk_targets_from_gamma(gamma, ctx))
  }
  targets <- lapply(ctx$comp$units, function(unit) {
    le_loo_target_from_gamma(gamma, ctx, kind = kind, unit = unit)
  })
  do.call(rbind, targets)
}

le_loo_target_from_gamma <- function(gamma, ctx,
                                     kind = c("jackknife", "ajk"), unit) {
  kind <- match.arg(kind)
  if (kind == "jackknife") {
    return(le_target_from_gamma(gamma, ctx, drop_unit = unit))
  }
  targets <- le_ajk_targets_from_gamma(gamma, ctx)
  targets[match(unit, ctx$comp$units), ]
}

le_ajk_targets_from_gamma <- function(gamma, ctx) {
  comp <- le_comp_with_gamma(ctx$comp, gamma)
  spec <- le_link_spec(comp, ctx)
  le_ajk_targets_analytic(
    delta_hat = spec$delta,
    Gi = spec$gmat_fn(spec$delta),
    A_full = spec$bread_fn(spec$delta),
    unit_bread_fn = function(unit) spec$unit_bread_fn(spec$delta, unit),
    target_fn = spec$target_fn,
    units = comp$units,
    unit_sets = as.list(seq_along(comp$units))
  )
}

le_target_from_gamma <- function(gamma, ctx, drop_unit = NULL) {
  comp <- le_comp_with_gamma(ctx$comp, gamma)
  if (!is.null(drop_unit)) {
    comp <- le_comp_drop_unit(comp, drop_unit)
  }
  spec <- le_link_spec(comp, ctx)
  spec$target_fn(spec$delta)
}

le_link_spec <- function(comp, ctx) {
  ref_idx <- ctx$ref_idx
  if (ctx$approach == "chain") {
    delta <- le_chain_delta(comp, type = ctx$method)
    return(list(
      delta = delta,
      target_fn = le_ref_target_fn(
        function(d) le_chain_target(d, type = ctx$method, Tn = comp$Tn),
        ref_idx),
      gmat_fn = function(d) le_chain_gmat(d, comp, type = ctx$method),
      bread_fn = function(d) le_chain_bread(d, comp, type = ctx$method),
      unit_bread_fn = function(d, unit) {
        le_chain_bread_unit(d, comp, type = ctx$method, unit = unit)
      }
    ))
  }
  omega <- le_phl_weights(comp, variant = ctx$phl_weights)
  delta <- le_phl_delta(comp, omega)
  list(
    delta = delta,
    target_fn = le_ref_target_fn(
      function(d) le_phl_target(d, Tn = comp$Tn), ref_idx),
    gmat_fn = function(d) le_phl_gmat(d, comp, omega),
    bread_fn = function(d) le_phl_bread(d, comp, omega),
    unit_bread_fn = function(d, unit) {
      le_phl_bread_unit(d, comp, omega, unit = unit)
    }
  )
}

le_gamma_vector <- function(comp) {
  x <- as.vector(t(as.matrix(comp$ip[c("a", "b")])))
  names(x) <- le_gamma_names(comp)
  x
}

le_gamma_names <- function(comp) {
  as.vector(t(cbind(
    paste0("g", comp$ip$group, ":", comp$ip$item, ":a"),
    paste0("g", comp$ip$group, ":", comp$ip$item, ":b")
  )))
}

le_comp_with_gamma <- function(comp, gamma) {
  gamma <- as.numeric(gamma)
  if (length(gamma) != 2L * nrow(comp$ip)) {
    stop("`vcov`/gamma length does not match link item parameters.",
         call. = FALSE)
  }
  ip <- comp$ip
  ip$a <- gamma[seq(1L, length(gamma), by = 2L)]
  ip$b <- gamma[seq(2L, length(gamma), by = 2L)]
  le_rebuild_comp(comp, ip = ip, units = comp$units)
}

le_comp_drop_unit <- function(comp, unit) {
  ip <- comp$ip[comp$ip$item != unit, , drop = FALSE]
  le_rebuild_comp(comp, ip = ip, units = setdiff(comp$units, unit))
}

le_rebuild_comp <- function(comp, ip, units) {
  by_group <- lapply(comp$groups, function(w) ip[ip$group == w, , drop = FALSE])
  if (any(vapply(by_group, nrow, integer(1)) == 0L)) {
    stop("Deleting a jackknife item removed all item parameters from a group.",
         call. = FALSE)
  }
  comp$ip <- ip
  comp$aL <- lapply(by_group, function(d) stats::setNames(d$a, d$item))
  comp$bL <- lapply(by_group, function(d) stats::setNames(d$b, d$item))
  comp$items_t <- lapply(by_group, function(d) d$item)
  comp$units <- units
  comp
}

le_vcov_matrix <- function(vcov, comp) {
  expected <- le_gamma_names(comp)
  p <- length(expected)
  if (is.matrix(vcov)) {
    V <- le_reorder_vcov(vcov, expected = expected,
                         label = "full `vcov`")
    return(V)
  }
  if (!is.list(vcov)) {
    stop("`vcov` must be a covariance matrix or a per-group list of ",
         "covariance matrices.", call. = FALSE)
  }
  V <- matrix(0, nrow = p, ncol = p, dimnames = list(expected, expected))
  for (ww in seq_along(comp$groups)) {
    group <- comp$groups[ww]
    Vw <- vcov[[as.character(group)]]
    if (is.null(Vw)) Vw <- vcov[[ww]]
    if (is.null(Vw)) {
      stop("Missing `vcov` matrix for group ", group, ".", call. = FALSE)
    }
    rows <- which(comp$ip$group == group)
    pos <- as.vector(t(cbind(2L * rows - 1L, 2L * rows)))
    items <- comp$ip$item[rows]
    expected_full <- expected[pos]
    expected_short <- as.vector(t(cbind(paste0(items, "_a"),
                                         paste0(items, "_b"))))
    V[pos, pos] <- le_reorder_vcov(Vw, expected = expected_full,
                                   fallback = expected_short,
                                   label = paste0("`vcov` for group ", group))
  }
  V
}

le_reorder_vcov <- function(V, expected, fallback = NULL, label = "`vcov`") {
  if (!is.matrix(V)) {
    stop(label, " must be a matrix.", call. = FALSE)
  }
  n <- length(expected)
  rn <- rownames(V)
  cn <- colnames(V)
  if (!is.null(rn) && !is.null(cn)) {
    if (all(expected %in% rn) && all(expected %in% cn)) {
      return(V[expected, expected, drop = FALSE])
    }
    if (!is.null(fallback) && all(fallback %in% rn) &&
        all(fallback %in% cn)) {
      out <- V[fallback, fallback, drop = FALSE]
      dimnames(out) <- list(expected, expected)
      return(out)
    }
  }
  if (!identical(dim(V), c(n, n))) {
    stop(label, " has dimension ", paste(dim(V), collapse = "x"),
         "; expected ", n, "x", n, ".", call. = FALSE)
  }
  dimnames(V) <- list(expected, expected)
  V
}

le_error_table <- function(V_SE, V_LE, V_LEbc, groups, method, n_reps,
                           ref_idx = NULL) {
  se <- le_sqrt_diag(V_SE)
  le <- le_sqrt_diag(V_LE)
  lebc <- le_sqrt_diag(V_LEbc)
  P <- length(groups) - 1L
  pick <- function(x) {
    out <- cbind(sigma = c(0, x[seq_len(P)]),
                 mu = c(0, x[P + seq_len(P)]))
    if (!is.null(ref_idx)) {
      out[1L, "sigma"] <- x[2L * P + 1L]
      out[1L, "mu"] <- x[2L * P + 2L]
      out[ref_idx, ] <- 0
    }
    out
  }
  se_m <- pick(se); le_m <- pick(le); lebc_m <- pick(lebc)
  te_m <- sqrt(se_m^2 + le_m^2)
  tebc_m <- sqrt(se_m^2 + lebc_m^2)
  out <- data.frame(
    group = groups,
    le_mu = le_m[, "mu"],
    le_sigma = le_m[, "sigma"],
    method = method,
    n_reps = n_reps,
    se_mu = se_m[, "mu"],
    se_sigma = se_m[, "sigma"],
    lebc_mu = lebc_m[, "mu"],
    lebc_sigma = lebc_m[, "sigma"],
    te_mu = te_m[, "mu"],
    te_sigma = te_m[, "sigma"],
    tebc_mu = tebc_m[, "mu"],
    tebc_sigma = tebc_m[, "sigma"]
  )
  rownames(out) <- NULL
  out
}

le_sqrt_diag <- function(V) {
  sqrt(pmax(diag(V), 0))
}

# expected sandwich, within-pattern item score covariances as meat
le_esw <- function(link) {
  if (identical(link$approach, "chain")) {
    if (!link$method %in% c("mgm", "mm")) {
      stop("Expected sandwich linking error for chain links is currently ",
           "only for method = \"mgm\" or \"mm\".", call. = FALSE)
    }
    return(le_esw_chain(link))
  }
  if (identical(link$approach, "joint")) {
    if (!identical(link$method, "phl")) {
      stop("Expected sandwich linking error for joint links is currently ",
           "only for method = \"phl\".", call. = FALSE)
    }
    return(le_esw_phl(link))
  }
  stop("Expected sandwich linking error is not currently available for ",
       "approach = \"", link$approach, "\".", call. = FALSE)
}

le_esw_chain <- function(link) {
  comp <- le_ipar_components(link, unit_mode = "chain")
  ref_idx <- le_ref_idx(link, comp$groups)
  delta <- le_chain_delta(comp, type = link$method)
  Gi <- le_chain_gmat(delta, comp, type = link$method)
  tgt <- le_chain_target(delta, type = link$method, Tn = comp$Tn)
  le_sandwich_from_mats(
    bread = le_chain_bread(delta, comp, type = link$method),
    B = le_expected_meat(Gi, le_unit_patterns(comp)),
    target_jac = le_rebase_jac(tgt, ref_idx) %*%
      le_chain_target_jac(delta, type = link$method, Tn = comp$Tn),
    groups = comp$groups,
    units = comp$units,
    method = "sandwich_esw",
    label = "Expected sandwich",
    ref_idx = ref_idx
  )
}

le_esw_phl <- function(link) {
  comp <- le_ipar_components(link, unit_mode = "joint")
  ref_idx <- le_ref_idx(link, comp$groups)
  omega <- le_phl_weights(comp, variant = le_phl_variant(link))
  delta <- le_phl_delta(comp, omega)
  Gi <- le_phl_gmat(delta, comp, omega)
  tgt <- le_phl_target(delta, Tn = comp$Tn)
  le_sandwich_from_mats(
    bread = le_phl_bread(delta, comp, omega),
    B = le_expected_meat(Gi, le_unit_patterns(comp)),
    target_jac = le_rebase_jac(tgt, ref_idx) %*%
      le_phl_target_jac(delta, Tn = comp$Tn),
    groups = comp$groups,
    units = comp$units,
    method = "sandwich_esw",
    label = "Expected sandwich",
    ref_idx = ref_idx
  )
}

# observed sandwich with analytic bread and target Jacobian
le_osw_analytic <- function(link, method = "sandwich_osw_analytic") {
  if (identical(link$approach, "chain")) {
    if (!link$method %in% c("mgm", "mm")) {
      stop("Analytic observed sandwich linking error for chain links is ",
           "currently available only for method = \"mgm\" or \"mm\".",
           call. = FALSE)
    }
    return(le_osw_analytic_chain(link, method = method))
  }
  if (identical(link$approach, "joint")) {
    if (!identical(link$method, "phl")) {
      stop("Analytic observed sandwich linking error for joint links is ",
           "currently available only for method = \"phl\".",
           call. = FALSE)
    }
    return(le_osw_analytic_phl(link, method = method))
  }
  stop("Analytic observed sandwich linking error is not currently available ",
       "for approach = \"", link$approach, "\".", call. = FALSE)
}

le_osw_analytic_chain <- function(link, method = "sandwich_osw_analytic") {
  comp <- le_ipar_components(link, unit_mode = "chain")
  ref_idx <- le_ref_idx(link, comp$groups)
  delta <- le_chain_delta(comp, type = link$method)
  Gi <- le_chain_gmat(delta, comp, type = link$method)
  tgt <- le_chain_target(delta, type = link$method, Tn = comp$Tn)
  le_sandwich_from_mats(
    bread = le_chain_bread(delta, comp, type = link$method),
    B = crossprod(Gi),
    target_jac = le_rebase_jac(tgt, ref_idx) %*%
      le_chain_target_jac(delta, type = link$method, Tn = comp$Tn),
    groups = comp$groups,
    units = comp$units,
    method = method,
    label = "Analytic observed sandwich",
    ref_idx = ref_idx
  )
}

le_osw_analytic_phl <- function(link, method = "sandwich_osw_analytic") {
  comp <- le_ipar_components(link, unit_mode = "joint")
  ref_idx <- le_ref_idx(link, comp$groups)
  omega <- le_phl_weights(comp, variant = le_phl_variant(link))
  delta <- le_phl_delta(comp, omega)
  Gi <- le_phl_gmat(delta, comp, omega)
  tgt <- le_phl_target(delta, Tn = comp$Tn)
  le_sandwich_from_mats(
    bread = le_phl_bread(delta, comp, omega),
    B = crossprod(Gi),
    target_jac = le_rebase_jac(tgt, ref_idx) %*%
      le_phl_target_jac(delta, Tn = comp$Tn),
    groups = comp$groups,
    units = comp$units,
    method = method,
    label = "Analytic observed sandwich",
    ref_idx = ref_idx
  )
}

le_sandwich_from_mats <- function(bread, B, target_jac, groups, units,
                                  method, label = "Sandwich",
                                  ref_idx = NULL) {
  Ainv <- tryCatch(
    solve(bread),
    error = function(e) {
      stop(label, " bread matrix is singular; linking error cannot be ",
           "estimated. Original error: ", conditionMessage(e),
           call. = FALSE)
    }
  )
  V <- Ainv %*% B %*% t(Ainv)
  Vt <- target_jac %*% V %*% t(target_jac)
  se <- sqrt(pmax(diag(Vt), 0))
  le_from_target_se(se, groups = groups, method = method,
                    n_reps = length(units), ref_idx = ref_idx)
}

# observed sandwich linking error
le_sandwich <- function(link, method) {
  correction <- switch(method,
    sandwich_osw = "osw",
    sandwich_bosw = "bosw"
  )
  if (identical(link$approach, "chain")) {
    if (!link$method %in% c("mgm", "mm")) {
      stop("Sandwich linking error for chain links is currently available ",
           "only for method = \"mgm\" or \"mm\".", call. = FALSE)
    }
    return(le_sandwich_chain(link, method, correction))
  }
  if (identical(link$approach, "joint")) {
    if (!identical(link$method, "phl")) {
      stop("Sandwich linking error for joint links is currently available ",
           "only for method = \"phl\".", call. = FALSE)
    }
    return(le_sandwich_phl(link, method, correction))
  }
  stop("Sandwich linking error is not currently available for approach = \"",
       link$approach, "\".", call. = FALSE)
}

le_sandwich_chain <- function(link, method, correction) {
  comp <- le_ipar_components(link, unit_mode = "chain")
  ref_idx <- le_ref_idx(link, comp$groups)
  delta <- le_chain_delta(comp, type = link$method)
  le_sandwich_core(
    delta_hat = delta,
    gmat_fn = function(d) le_chain_gmat(d, comp, type = link$method),
    target_fn = le_ref_target_fn(
      function(d) le_chain_target(d, type = link$method, Tn = comp$Tn),
      ref_idx),
    groups = comp$groups,
    units = comp$units,
    method = method,
    correction = correction,
    ref_idx = ref_idx
  )
}

le_sandwich_phl <- function(link, method, correction) {
  comp <- le_ipar_components(link, unit_mode = "joint")
  ref_idx <- le_ref_idx(link, comp$groups)
  omega <- le_phl_weights(comp, variant = le_phl_variant(link))
  delta <- le_phl_delta(comp, omega)
  le_sandwich_core(
    delta_hat = delta,
    gmat_fn = function(d) le_phl_gmat(d, comp, omega),
    target_fn = le_ref_target_fn(
      function(d) le_phl_target(d, Tn = comp$Tn), ref_idx),
    groups = comp$groups,
    units = comp$units,
    method = method,
    correction = correction,
    ref_idx = ref_idx
  )
}

# bias-corrected sandwich linking errors, bias-corrected meat
# B - H V_diag H^T (Robitzsch, 2024, Stats, Eqs. 35-36)
le_sandwich_bc <- function(link, vcov, method) {
  if (is.null(vcov)) {
    stop("Bias-corrected sandwich linking error needs `vcov`, the ",
         "item-parameter covariance matrix/list.", call. = FALSE)
  }
  ctx <- le_bc_context(link, what = "Bias-corrected sandwich linking error")
  comp <- ctx$comp
  full_ip_key <- paste(comp$ip$group, comp$ip$item)
  n_gamma <- 2L * nrow(comp$ip)
  Vg <- le_vcov_matrix(vcov, comp)
  spec <- le_link_spec(comp, ctx)
  delta <- spec$delta
  if (ctx$approach == "chain") {
    A <- le_chain_bread(delta, comp, type = ctx$method)
    H <- le_chain_hgamma(delta, comp, type = ctx$method,
                         full_ip_key = full_ip_key, n_gamma = n_gamma)
    Jt <- le_rebase_jac(le_chain_target(delta, type = ctx$method,
                                        Tn = comp$Tn), ctx$ref_idx) %*%
      le_chain_target_jac(delta, type = ctx$method, Tn = comp$Tn)
    Gi <- le_chain_gmat(delta, comp, type = ctx$method)
  } else {
    omega <- le_phl_weights(comp, variant = ctx$phl_weights)
    A <- le_phl_bread(delta, comp, omega)
    H <- le_phl_hgamma(delta, comp, omega,
                       full_ip_key = full_ip_key, n_gamma = n_gamma)
    Jt <- le_rebase_jac(le_phl_target(delta, Tn = comp$Tn),
                        ctx$ref_idx) %*%
      le_phl_target_jac(delta, Tn = comp$Tn)
    Gi <- le_phl_gmat(delta, comp, omega)
  }
  Ainv <- tryCatch(
    solve(A),
    error = function(e) {
      stop("Sandwich bread matrix is singular; linking error cannot be ",
           "estimated. Original error: ", conditionMessage(e),
           call. = FALSE)
    }
  )
  B <- if (identical(method, "sandwich_esw_bc")) {
    le_expected_meat(Gi, le_unit_patterns(comp))
  } else {
    crossprod(Gi)
  }
  n_units <- length(comp$units)
  fac <- if (identical(method, "sandwich_bosw_bc")) {
    n_units / (n_units - 1)
  } else {
    1
  }
  U <- -Jt %*% Ainv %*% H
  V_SE <- U %*% Vg %*% t(U)
  Vg_diag <- matrix(0, nrow(Vg), ncol(Vg))
  for (it in unique(comp$ip$item)) {
    rows <- which(comp$ip$item == it)
    idx <- as.integer(rbind(2L * rows - 1L, 2L * rows))
    Vg_diag[idx, idx] <- Vg[idx, idx]
  }
  D_tilde <- H %*% Vg_diag %*% t(H)
  core <- Jt %*% Ainv
  V_LE <- fac * core %*% B %*% t(core)
  V_LEbc <- fac * core %*% (B - D_tilde) %*% t(core)
  le_error_table(V_SE = V_SE, V_LE = V_LE, V_LEbc = V_LEbc,
                 groups = comp$groups, method = method,
                 n_reps = n_units, ref_idx = ctx$ref_idx)
}

le_ipar_components <- function(link, unit_mode = c("chain", "joint")) {
  unit_mode <- match.arg(unit_mode)
  ip <- as_ipars(link$ipars)
  if (!is.null(link$anchor)) {
    ip <- ip[ip$item %in% link$anchor, , drop = FALSE]
  }
  groups <- link$trend$group
  by_group <- lapply(groups, function(w) ip[ip$group == w, , drop = FALSE])
  if (any(vapply(by_group, nrow, integer(1)) == 0L)) {
    stop("Sandwich linking error needs item parameters for every linked group.",
         call. = FALSE)
  }
  aL <- lapply(by_group, function(d) stats::setNames(d$a, d$item))
  bL <- lapply(by_group, function(d) stats::setNames(d$b, d$item))
  items_t <- lapply(by_group, function(d) d$item)
  if (any(unlist(aL, use.names = FALSE) <= 0)) {
    stop("Non-positive discrimination(s) in item parameters; ",
         "remove degenerate items before estimating linking error.",
         call. = FALSE)
  }

  units <- if (unit_mode == "chain") {
    Reduce(union, lapply(seq_len(length(groups) - 1L), function(t) {
      intersect(items_t[[t]], items_t[[t + 1L]])
    }))
  } else {
    presence <- table(ip$item)
    names(presence)[presence >= 2L]
  }
  units <- sort(unique(units))
  if (length(units) < 2L) {
    stop("Sandwich linking error needs at least two link items.",
         call. = FALSE)
  }
  list(ip = ip, groups = groups, Tn = length(groups), aL = aL, bL = bL,
       items_t = items_t, units = units)
}

le_unit_patterns <- function(comp) {
  vapply(comp$units, function(item) {
    paste(as.integer(vapply(comp$items_t, function(items) item %in% items,
                            logical(1))), collapse = "")
  }, character(1))
}

le_expected_meat <- function(Gi, patterns) {
  Gi <- as.matrix(Gi)
  if (nrow(Gi) != length(patterns)) {
    stop("Expected sandwich meat needs one item pattern per score row.",
         call. = FALSE)
  }
  B <- matrix(0, ncol = ncol(Gi), nrow = ncol(Gi))
  for (pat in unique(patterns)) {
    idx <- which(patterns == pat)
    if (length(idx) < 2L) {
      stop("Expected sandwich linking error needs at least two link items ",
           "per item-presence pattern. Pattern ", pat, " has only one.",
           call. = FALSE)
    }
    Gp <- Gi[idx, , drop = FALSE]
    C <- stats::cov(Gp)
    if (!is.matrix(C)) C <- matrix(C, nrow = 1L, ncol = 1L)
    C[!is.finite(C)] <- 0
    B <- B + length(idx) * C
  }
  B
}

le_sandwich_core <- function(delta_hat, gmat_fn, target_fn, groups, units,
                             method, correction, ref_idx = NULL) {
  Gsum <- function(d) colSums(gmat_fn(d))
  A <- le_numjac(Gsum, delta_hat)
  Gi <- gmat_fn(delta_hat)
  B <- crossprod(Gi)
  Ainv <- tryCatch(
    solve(A),
    error = function(e) {
      stop("Sandwich bread matrix is singular; linking error cannot be ",
           "estimated. Original error: ", conditionMessage(e),
           call. = FALSE)
    }
  )
  V <- Ainv %*% B %*% t(Ainv)
  Jt <- le_numjac(target_fn, delta_hat)
  Vt <- Jt %*% V %*% t(Jt)
  se <- sqrt(pmax(diag(Vt), 0))
  if (correction == "bosw") {
    se <- sqrt(length(units) / (length(units) - 1)) * se
  }
  le_from_target_se(se, groups = groups, method = method,
                    n_reps = length(units), ref_idx = ref_idx)
}

le_from_target_se <- function(se, groups, method, n_reps, ref_idx = NULL) {
  P <- length(groups) - 1L
  le <- data.frame(
    group = groups,
    le_mu = c(0, se[P + seq_len(P)]),
    le_sigma = c(0, se[seq_len(P)]),
    method = method,
    n_reps = n_reps
  )
  if (!is.null(ref_idx)) {
    # rebased targets, group 1 takes the appended slots and the
    # reference group is error free by construction
    le$le_sigma[1L] <- se[2L * P + 1L]
    le$le_mu[1L] <- se[2L * P + 2L]
    le$le_sigma[ref_idx] <- 0
    le$le_mu[ref_idx] <- 0
  }
  rownames(le) <- NULL
  le
}

# reference group support, rebasing only changes the projection
# side, the rebased vector appends the group-1 targets

# position of the reference group in the groups vector, NULL = default
le_ref_idx <- function(link, groups) {
  if (is.null(link$ref)) return(NULL)
  r <- match(link$ref, groups)
  if (is.na(r)) {
    stop("`ref` must be one of the linked groups (",
         paste(groups, collapse = ", "), "); got ", link$ref, ".",
         call. = FALSE)
  }
  if (r == 1L) return(NULL)
  r
}

le_rebase_target_vec <- function(target, ref_idx) {
  if (is.null(ref_idx)) return(target)
  P <- length(target) / 2L
  sigma <- c(1, target[seq_len(P)])
  mu <- c(0, target[P + seq_len(P)])
  s_r <- sigma[ref_idx]
  m_r <- mu[ref_idx]
  c(sigma[-1L] / s_r, (mu[-1L] - m_r) / s_r, 1 / s_r, -m_r / s_r)
}

le_ref_target_fn <- function(target_fn, ref_idx) {
  if (is.null(ref_idx)) return(target_fn)
  function(d) le_rebase_target_vec(target_fn(d), ref_idx)
}

# Jacobian of the rebasing map with respect to the original targets
# (2P + 2 rows, 2P columns), zero rows for t = r
le_rebase_jac <- function(target, ref_idx) {
  P <- length(target) / 2L
  if (is.null(ref_idx)) return(diag(2L * P))
  sigma <- c(1, target[seq_len(P)])
  mu <- c(0, target[P + seq_len(P)])
  s_r <- sigma[ref_idx]
  m_r <- mu[ref_idx]
  cs <- ref_idx - 1L        # Spalte von sigma_r
  cm <- P + ref_idx - 1L    # Spalte von mu_r
  R <- matrix(0, 2L * P + 2L, 2L * P)
  for (t in seq_len(P)) {
    R[t, t] <- 1 / s_r
    R[t, cs] <- R[t, cs] - sigma[t + 1L] / s_r^2
    R[P + t, P + t] <- 1 / s_r
    R[P + t, cm] <- R[P + t, cm] - 1 / s_r
    R[P + t, cs] <- R[P + t, cs] - (mu[t + 1L] - m_r) / s_r^2
  }
  R[2L * P + 1L, cs] <- -1 / s_r^2
  R[2L * P + 2L, cm] <- -1 / s_r
  R[2L * P + 2L, cs] <- m_r / s_r^2
  R
}

le_numjac <- function(f, x, h = 1e-5) {
  y0 <- f(x)
  J <- matrix(NA_real_, nrow = length(y0), ncol = length(x))
  for (j in seq_along(x)) {
    step <- h * (1 + abs(x[j]))
    xp <- xm <- x
    xp[j] <- xp[j] + step
    xm[j] <- xm[j] - step
    J[, j] <- (f(xp) - f(xm)) / (2 * step)
  }
  J
}

le_chain_delta <- function(comp, type = c("mgm", "mm")) {
  type <- match.arg(type)
  P <- comp$Tn - 1L
  delta <- numeric(2L * P)
  for (step in seq_len(P)) {
    com <- intersect(comp$items_t[[step]], comp$items_t[[step + 1L]])
    a0 <- comp$aL[[step]][com]
    a1 <- comp$aL[[step + 1L]][com]
    b0 <- comp$bL[[step]][com]
    b1 <- comp$bL[[step + 1L]][com]
    ratio <- if (type == "mgm") {
      exp(mean(log(a1) - log(a0)))
    } else {
      sum(a1) / sum(a0)
    }
    delta[2L * step - 1L] <- if (type == "mgm") log(ratio) else ratio
    delta[2L * step] <- mean(b0 - ratio * b1)
  }
  delta
}

le_chain_target <- function(delta, type = c("mgm", "mm"), Tn) {
  type <- match.arg(type)
  P <- Tn - 1L
  scale_idx <- 2L * seq_len(P) - 1L
  mean_idx <- 2L * seq_len(P)
  sigma_step <- if (type == "mgm") exp(delta[scale_idx]) else delta[scale_idx]
  mu_step <- delta[mean_idx]
  sigma <- numeric(Tn)
  mu <- numeric(Tn)
  sigma[1L] <- 1
  for (t in 2:Tn) {
    sigma[t] <- sigma[t - 1L] * sigma_step[t - 1L]
    mu[t] <- mu[t - 1L] + sigma[t - 1L] * mu_step[t - 1L]
  }
  c(sigma[-1L], mu[-1L])
}

le_chain_target_jac <- function(delta, type = c("mgm", "mm"), Tn) {
  type <- match.arg(type)
  P <- Tn - 1L
  scale_idx <- 2L * seq_len(P) - 1L
  mean_idx <- 2L * seq_len(P)
  sigma_step <- if (type == "mgm") exp(delta[scale_idx]) else delta[scale_idx]
  dlog_step <- if (type == "mgm") rep(1, P) else 1 / delta[scale_idx]
  mu_step <- delta[mean_idx]

  sigma <- numeric(Tn)
  mu <- numeric(Tn)
  sigma[1L] <- 1
  for (t in 2:Tn) {
    sigma[t] <- sigma[t - 1L] * sigma_step[t - 1L]
    mu[t] <- mu[t - 1L] + sigma[t - 1L] * mu_step[t - 1L]
  }

  J <- matrix(0, nrow = 2L * P, ncol = 2L * P)
  for (out in seq_len(P)) {
    group <- out + 1L
    for (step in seq_len(out)) {
      J[out, scale_idx[step]] <- sigma[group] * dlog_step[step]
    }
  }
  for (out in seq_len(P)) {
    group <- out + 1L
    row <- P + out
    for (step in seq_len(out)) {
      J[row, mean_idx[step]] <- sigma[step]
      if (step < out) {
        contrib <- sum(sigma[(step + 1L):out] * mu_step[(step + 1L):out])
        J[row, scale_idx[step]] <- dlog_step[step] * contrib
      }
    }
  }
  J
}

le_chain_bread <- function(delta, comp, type = c("mgm", "mm")) {
  type <- match.arg(type)
  P <- comp$Tn - 1L
  A <- matrix(0, nrow = 2L * P, ncol = 2L * P)
  for (step in seq_len(P)) {
    sp <- 2L * step - 1L
    mp <- sp + 1L
    com <- intersect(comp$items_t[[step]], comp$items_t[[step + 1L]])
    com <- intersect(com, comp$units)
    a0 <- comp$aL[[step]][com]
    b1 <- comp$bL[[step + 1L]][com]
    if (type == "mgm") {
      A[sp, sp] <- -length(com)
      A[mp, sp] <- sum(exp(delta[sp]) * b1)
      A[mp, mp] <- length(com)
    } else {
      A[sp, sp] <- sum(a0)
      A[mp, sp] <- sum(b1)
      A[mp, mp] <- length(com)
    }
  }
  A
}

le_chain_bread_unit <- function(delta, comp, type = c("mgm", "mm"), unit) {
  type <- match.arg(type)
  P <- comp$Tn - 1L
  A <- matrix(0, nrow = 2L * P, ncol = 2L * P)
  for (step in seq_len(P)) {
    com <- intersect(comp$items_t[[step]], comp$items_t[[step + 1L]])
    if (!unit %in% com || !unit %in% comp$units) {
      next
    }
    sp <- 2L * step - 1L
    mp <- sp + 1L
    a0 <- comp$aL[[step]][[unit]]
    b1 <- comp$bL[[step + 1L]][[unit]]
    if (type == "mgm") {
      A[sp, sp] <- -1
      A[mp, sp] <- exp(delta[sp]) * b1
      A[mp, mp] <- 1
    } else {
      A[sp, sp] <- a0
      A[mp, sp] <- b1
      A[mp, mp] <- 1
    }
  }
  A
}

le_chain_gmat <- function(delta, comp, type = c("mgm", "mm")) {
  type <- match.arg(type)
  P <- comp$Tn - 1L
  M <- matrix(0, nrow = length(comp$units), ncol = 2L * P,
              dimnames = list(comp$units, NULL))
  for (step in seq_len(P)) {
    sp <- 2L * step - 1L
    mp <- sp + 1L
    sc <- delta[sp]
    mn <- delta[mp]
    com <- intersect(comp$items_t[[step]], comp$items_t[[step + 1L]])
    com <- intersect(com, comp$units)
    ix <- match(com, comp$units)
    a0 <- comp$aL[[step]][com]
    a1 <- comp$aL[[step + 1L]][com]
    b0 <- comp$bL[[step]][com]
    b1 <- comp$bL[[step + 1L]][com]
    if (type == "mgm") {
      M[ix, sp] <- log(a1) - log(a0) - sc
      M[ix, mp] <- exp(sc) * b1 - b0 + mn
    } else {
      M[ix, sp] <- sc * a0 - a1
      M[ix, mp] <- sc * b1 - b0 + mn
    }
  }
  M
}

# item_weights setting of a joint PHL link
le_phl_variant <- function(link) {
  # exact [[, `$` could partial-match another control entry
  w <- link$control[["item_weights"]]
  if (is.null(w)) "inverse_admin" else w
}

le_phl_weights <- function(comp, variant = c("inverse_admin", "uniform")) {
  variant <- match.arg(variant)
  presence <- table(unlist(comp$items_t, use.names = FALSE))
  if (variant == "uniform") {
    stats::setNames(rep(1, length(presence)), names(presence))
  } else {
    stats::setNames(as.numeric(1 / presence), names(presence))
  }
}

le_phl_delta <- function(comp, omega) {
  P <- comp$Tn - 1L
  idx <- function(t) t - 1L
  XtWX_s <- matrix(0, P, P)
  XtWy_s <- numeric(P)
  rows <- list()
  r <- 0L
  for (g in seq_len(comp$Tn - 1L)) {
    for (h in (g + 1L):comp$Tn) {
      com <- intersect(comp$items_t[[g]], comp$items_t[[h]])
      for (it in com) {
        w <- omega[[it]]
        y <- log(comp$aL[[g]][[it]]) - log(comp$aL[[h]][[it]])
        x <- numeric(P)
        if (g >= 2L) x[idx(g)] <- x[idx(g)] + 1
        if (h >= 2L) x[idx(h)] <- x[idx(h)] - 1
        XtWX_s <- XtWX_s + w * (x %o% x)
        XtWy_s <- XtWy_s + w * x * y
        r <- r + 1L
        rows[[r]] <- list(g = g, h = h, it = it, w = w, x = x)
      }
    }
  }
  if (r == 0L) {
    stop("No common items between any pair of groups; PHL is not possible.",
         call. = FALSE)
  }
  s_free <- tryCatch(
    solve(XtWX_s, XtWy_s),
    error = function(e) {
      stop("PHL scale system is singular; the group link graph may be ",
           "disconnected. Original error: ", conditionMessage(e),
           call. = FALSE)
    }
  )
  sigma <- exp(c(0, s_free))
  XtWy_m <- numeric(P)
  for (rr in rows) {
    y <- -(sigma[rr$g] * comp$bL[[rr$g]][[rr$it]] -
           sigma[rr$h] * comp$bL[[rr$h]][[rr$it]])
    XtWy_m <- XtWy_m + rr$w * rr$x * y
  }
  mu_free <- tryCatch(
    solve(XtWX_s, XtWy_m),
    error = function(e) {
      stop("PHL mean system is singular; the group link graph may be ",
           "disconnected. Original error: ", conditionMessage(e),
           call. = FALSE)
    }
  )
  c(s_free, mu_free)
}

le_phl_target <- function(delta, Tn) {
  P <- Tn - 1L
  c(exp(delta[seq_len(P)]), delta[P + seq_len(P)])
}

le_phl_target_jac <- function(delta, Tn) {
  P <- Tn - 1L
  J <- matrix(0, nrow = 2L * P, ncol = 2L * P)
  J[cbind(seq_len(P), seq_len(P))] <- exp(delta[seq_len(P)])
  J[cbind(P + seq_len(P), P + seq_len(P))] <- 1
  J
}

le_phl_bread <- function(delta, comp, omega) {
  P <- comp$Tn - 1L
  s <- c(0, delta[seq_len(P)])
  sigma <- exp(s)
  A <- matrix(0, nrow = 2L * P, ncol = 2L * P)
  for (ii in seq_along(comp$units)) {
    it <- comp$units[ii]
    tps <- which(vapply(comp$items_t, function(v) it %in% v, logical(1)))
    if (length(tps) < 2L) next
    w <- omega[[it]]
    for (k in tps[tps >= 2L]) {
      oth <- setdiff(tps, k)
      r_s <- k - 1L
      r_m <- P + k - 1L
      A[r_s, k - 1L] <- A[r_s, k - 1L] + w * length(oth)
      A[r_m, k - 1L] <- A[r_m, k - 1L] +
        w * length(oth) * sigma[k] * comp$bL[[k]][[it]]
      A[r_m, P + k - 1L] <- A[r_m, P + k - 1L] + w * length(oth)
      for (h in oth[oth >= 2L]) {
        A[r_s, h - 1L] <- A[r_s, h - 1L] - w
        A[r_m, h - 1L] <- A[r_m, h - 1L] -
          w * sigma[h] * comp$bL[[h]][[it]]
        A[r_m, P + h - 1L] <- A[r_m, P + h - 1L] - w
      }
    }
  }
  A
}

le_phl_bread_unit <- function(delta, comp, omega, unit) {
  P <- comp$Tn - 1L
  s <- c(0, delta[seq_len(P)])
  sigma <- exp(s)
  A <- matrix(0, nrow = 2L * P, ncol = 2L * P)
  if (!unit %in% comp$units) {
    return(A)
  }
  tps <- which(vapply(comp$items_t, function(v) unit %in% v, logical(1)))
  if (length(tps) < 2L) {
    return(A)
  }
  w <- omega[[unit]]
  for (k in tps[tps >= 2L]) {
    oth <- setdiff(tps, k)
    r_s <- k - 1L
    r_m <- P + k - 1L
    A[r_s, k - 1L] <- A[r_s, k - 1L] + w * length(oth)
    A[r_m, k - 1L] <- A[r_m, k - 1L] +
      w * length(oth) * sigma[k] * comp$bL[[k]][[unit]]
    A[r_m, P + k - 1L] <- A[r_m, P + k - 1L] + w * length(oth)
    for (h in oth[oth >= 2L]) {
      A[r_s, h - 1L] <- A[r_s, h - 1L] - w
      A[r_m, h - 1L] <- A[r_m, h - 1L] -
        w * sigma[h] * comp$bL[[h]][[unit]]
      A[r_m, P + h - 1L] <- A[r_m, P + h - 1L] - w
    }
  }
  A
}

le_phl_gmat <- function(delta, comp, omega) {
  P <- comp$Tn - 1L
  s <- c(0, delta[seq_len(P)])
  mu <- c(0, delta[P + seq_len(P)])
  sigma <- exp(s)
  M <- matrix(0, nrow = length(comp$units), ncol = 2L * P,
              dimnames = list(comp$units, NULL))
  for (ii in seq_along(comp$units)) {
    it <- comp$units[ii]
    tps <- which(vapply(comp$items_t, function(v) it %in% v, logical(1)))
    if (length(tps) < 2L) next
    w <- omega[[it]]
    for (k in tps[tps >= 2L]) {
      oth <- setdiff(tps, k)
      a_oth <- vapply(oth, function(h) log(comp$aL[[h]][[it]]), numeric(1))
      b_oth <- vapply(oth, function(h) sigma[h] * comp$bL[[h]][[it]],
                       numeric(1))
      M[ii, k - 1L] <- w * sum((s[k] - s[oth]) -
                                 (log(comp$aL[[k]][[it]]) - a_oth))
      M[ii, P + k - 1L] <- w * sum(
        (sigma[k] * comp$bL[[k]][[it]] - b_oth) + (mu[k] - mu[oth])
      )
    }
  }
  M
}
