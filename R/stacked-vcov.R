# R/stacked-vcov.R
# joint sandwich covariance of separately calibrated groups for
# dependent person samples

# shared persons per group pair
vcov_overlap_matrix <- function(ids) {
  Tn <- length(ids)
  overlap <- matrix(0L, Tn, Tn,
                    dimnames = list(paste0("g", seq_len(Tn)),
                                    paste0("g", seq_len(Tn))))
  for (t in seq_len(Tn)) {
    for (u in seq_len(Tn)) {
      overlap[t, u] <- length(intersect(ids[[t]], ids[[u]]))
    }
  }
  overlap
}

# per-person score matrix of the marginal loglikelihood via the
# Fisher identity
em_person_scores <- function(fit) {
  em_vcov_supported(fit)
  ed <- em_estep_data(fit$dat, rep(1, nrow(fit$dat)))
  P <- em_irf_matrix(fit$a, fit$nu, fit$theta)
  es <- em_estep(ed, log(P), log(1 - P), log(fit$pi),
                 return_posterior = TRUE)
  W <- es$posterior                       # N x K, rows sum to 1
  WP <- W %*% t(P)                        # N x J: E[P_j(theta) | x_n]
  S_v <- ed$Obs * WP - ed$Y1
  if (fit$model != "2PL") {
    colnames(S_v) <- fit$item
    return(S_v)
  }
  Wt <- as.vector(W %*% fit$theta)        # N: E[theta | x_n]
  WPt <- W %*% t(sweep(P, 2L, fit$theta, "*"))
  S_a <- ed$Y1 * Wt - ed$Obs * WPt
  J <- length(fit$item)
  S <- matrix(0, nrow(S_a), 2L * J)
  S[, seq(1L, 2L * J, 2L)] <- S_a
  S[, seq(2L, 2L * J, 2L)] <- S_v
  colnames(S) <- as.vector(t(cbind(paste0(fit$item, "_a"),
                                   paste0(fit$item, "_v"))))
  S
}

stacked_vcov <- function(cal, ids, h = 1e-5) {
  if (!inherits(cal, "irtlink_calib"))
    stop("`cal` must be an irtlink_calib object.", call. = FALSE)
  if (!identical(cal$calibration, "separate"))
    stop("stacked_vcov requires calibration = \"separate\".", call. = FALSE)
  if (!identical(cal$engine, "em"))
    stop("stacked_vcov requires engine = \"em\".", call. = FALSE)
  fits <- cal$models
  if (is.null(fits))
    stop("stacked_vcov needs the stored EM fits; call calibrate() with ",
         "keep_models = TRUE.", call. = FALSE)
  Tn <- length(fits)
  if (!is.list(ids) || length(ids) != Tn)
    stop("`ids` must be a list with one identifier vector per group.",
         call. = FALSE)
  ids <- lapply(ids, as.character)
  for (t in seq_len(Tn)) {
    if (length(ids[[t]]) != nrow(fits[[t]]$dat))
      stop("`ids[[", t, "]]` must have one identifier per row of the ",
           "group-", t, " response data (", nrow(fits[[t]]$dat),
           " rows).", call. = FALSE)
    if (anyDuplicated(ids[[t]]))
      stop("`ids[[", t, "]]` contains duplicated identifiers.",
           call. = FALSE)
  }

  # per-group pieces, person scores, observed information, delta map
  S_list <- lapply(fits, em_person_scores)
  I_list <- lapply(fits, em_obs_info_ipars, h = h)
  D_list <- lapply(fits, em_delta_ab)
  p_vec <- vapply(I_list, nrow, integer(1))
  offs <- c(0L, cumsum(p_vec))
  p_tot <- offs[Tn + 1L]
  ab_names <- unlist(lapply(seq_len(Tn), function(t) {
    as.vector(t(cbind(paste0("g", t, ":", fits[[t]]$item, ":a"),
                      paste0("g", t, ":", fits[[t]]$item, ":b"))))
  }))

  # stacked score matrix, zero rows where a person missed the group
  all_ids <- unique(unlist(ids))
  S_full <- matrix(0, length(all_ids), p_tot)
  for (t in seq_len(Tn)) {
    rows <- match(ids[[t]], all_ids)
    S_full[rows, offs[t] + seq_len(p_vec[t])] <- S_list[[t]]
  }
  M <- crossprod(S_full)

  # bread inverse and delta map, block-diagonal per group
  # (per-group pinv reproduces em_vcov_ipars on the diagonal blocks)
  Kinv <- matrix(0, p_tot, p_tot)
  D_full <- matrix(0, 2L * sum(vapply(fits, function(f) length(f$item),
                                      integer(1))), p_tot)
  r0 <- 0L
  for (t in seq_len(Tn)) {
    idx <- offs[t] + seq_len(p_vec[t])
    Kinv[idx, idx] <- pinv(I_list[[t]])
    rt <- nrow(D_list[[t]])
    D_full[r0 + seq_len(rt), idx] <- D_list[[t]]
    r0 <- r0 + rt
  }

  to_ab <- function(V_coef) {
    V <- D_full %*% V_coef %*% t(D_full)
    V <- (V + t(V)) / 2
    dimnames(V) <- list(ab_names, ab_names)
    V
  }
  M_nocross <- matrix(0, p_tot, p_tot)
  for (t in seq_len(Tn)) {
    idx <- offs[t] + seq_len(p_vec[t])
    M_nocross[idx, idx] <- M[idx, idx]
  }
  vcov_dep <- to_ab(Kinv %*% M %*% Kinv)
  vcov_nocross <- to_ab(Kinv %*% M_nocross %*% Kinv)
  vcov_ind <- to_ab(Kinv)

  # per-item blocks across the groups administering the item
  items_all <- sort(unique(unlist(lapply(fits, `[[`, "item"))))
  vcov_items <- lapply(items_all, function(it) {
    nm <- ab_names[grepl(paste0(":", it, ":"), ab_names, fixed = TRUE)]
    vcov_dep[nm, nm, drop = FALSE]
  })
  names(vcov_items) <- items_all

  overlap <- vcov_overlap_matrix(ids)
  score_check <- vapply(S_list, function(S) max(abs(colSums(S))),
                        numeric(1))
  names(score_check) <- paste0("g", seq_len(Tn))

  structure(
    list(vcov = vcov_dep, vcov_nocross = vcov_nocross,
         vcov_independence = vcov_ind, vcov_items = vcov_items,
         overlap = overlap, score_check = score_check,
         n_groups = Tn, call = match.call()),
    class = "irtlink_stacked_vcov"
  )
}

print.irtlink_stacked_vcov <- function(x, ...) {
  cat("Stacked sandwich covariance of the item parameters (",
      x$n_groups, " groups, ", nrow(x$vcov), " parameters)\n", sep = "")
  cat("Shared persons per group pair:\n")
  print(x$overlap)
  cat("Max |score column sum| per group (gradient check):\n")
  print(signif(x$score_check, 3))
  invisible(x)
}
