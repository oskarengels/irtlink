# R/sim_trend_data.R
BASE_ITEMS <- data.frame(
  a = c(0.73, 1.25, 1.20, 1.47, 0.97, 1.38, 1.05, 1.14, 1.15, 0.67),
  b = c(-1.31, 1.44, -1.20, 0.10, 0.10, -0.74, 1.48, -0.61, 0.82, -0.07)
)

sim_trend_data <- function(n_groups, N, I = 20,
                           mu = NULL, sigma = NULL,
                           overlap = 0.4,
                           design = c("all", "successive"),
                           dif = c("none", "balanced", "unbalanced", "random"),
                           dif_pct = 0.3, dif_effect = 0.5, dif_sd = 0.3,
                           dif_par = c("b", "a", "both"),
                           dif_effect_a = 0.25, dif_sd_a = 0.15,
                           person_overlap = 0, person_cor = 0.7,
                           model = c("2PL", "1PL", "GPCM"),
                           n_cat = 3) {
  design <- match.arg(design)
  dif <- match.arg(dif)
  dif_par <- match.arg(dif_par)
  model <- match.arg(model)
  if (dif_par != "b" && model == "1PL") {
    stop("dif_par = \"", dif_par, "\" needs a model with free ",
         "discriminations (\"2PL\" or \"GPCM\").", call. = FALSE)
  }
  if (model == "GPCM") {
    if (length(n_cat) != 1L || n_cat != round(n_cat) || n_cat < 2 ||
        n_cat > 9) {
      stop("`n_cat` must be a single integer between 2 and 9.",
           call. = FALSE)
    }
  }
  n_thr <- if (model == "GPCM") as.integer(n_cat) - 1L else 1L
  if (person_overlap < 0 || person_overlap > 1) {
    stop("`person_overlap` must be in [0, 1].", call. = FALSE)
  }
  if (abs(person_cor) >= 1) {
    stop("`person_cor` must lie strictly between -1 and 1.", call. = FALSE)
  }
  if (n_groups < 2) {
    stop("`n_groups` must specify at least two groups.", call. = FALSE)
  }
  if (is.null(mu)) mu <- 0.3 * (seq_len(n_groups) - 1)
  if (is.null(sigma)) sigma <- rep(1, n_groups)
  if (length(mu) != n_groups) {
    stop("`mu` must have length `n_groups` (", n_groups, "); got ",
         length(mu), ".", call. = FALSE)
  }
  if (length(sigma) != n_groups) {
    stop("`sigma` must have length `n_groups` (", n_groups, "); got ",
         length(sigma), ".", call. = FALSE)
  }
  if (any(sigma <= 0)) {
    stop("`sigma` values must be positive.", call. = FALSE)
  }

  if (overlap <= 0 || overlap > 1) {
    stop("`overlap` must be in (0, 1].", call. = FALSE)
  }

  item_sets <- build_design_items(n_groups, I, overlap, design)
  pool <- sort(unique(unlist(item_sets)))
  idx <- ((seq_along(pool) - 1L) %% nrow(BASE_ITEMS)) + 1L
  true_ipars <- data.frame(
    item = pool,
    a = if (model == "1PL") rep(1, length(pool)) else BASE_ITEMS$a[idx],
    b = BASE_ITEMS$b[idx]
  )
  # GPCM thresholds b + tau_v with centered, equally spaced tau
  tau <- NULL
  if (model == "GPCM") {
    tau_v <- if (n_thr == 1L) 0 else seq(-0.6, 0.6, length.out = n_thr)
    tau <- matrix(tau_v, nrow = length(pool), ncol = n_thr, byrow = TRUE)
    colnames(tau) <- paste0("tau", seq_len(n_thr))
    true_ipars <- cbind(true_ipars, as.data.frame(tau))
  }

  # DIF shifts on link items at t >= 2, dif_par picks difficulties,
  # log discriminations, or both
  linkable <- names(which(table(unlist(item_sets)) >= 2))
  dif_items <- character(0)
  drift <- matrix(0, nrow = length(pool), ncol = n_groups,
                  dimnames = list(pool, NULL))
  drift_a <- drift
  if (dif != "none") {
    n_aff <- round(length(linkable) * dif_pct)
    if (n_aff < 1) {
      stop("`dif_pct` must select at least one affected item ",
           "(round(n_linkable * dif_pct) >= 1).", call. = FALSE)
    }
    dif_items <- sort(linkable)[seq_len(n_aff)]
    if (dif == "balanced" && n_aff %% 2 != 0) {
      warning("balanced DIF with an odd number of affected items (", n_aff,
              ") leaves a small net mean shift; use an even count for a ",
              "fully symmetric design.", call. = FALSE)
    }
    for (t in 2:n_groups) {
      for (k in seq_along(dif_items)) {
        it <- dif_items[k]
        if (dif_par %in% c("b", "both")) {
          drift[it, t] <- switch(dif,
            unbalanced = dif_effect,
            balanced = if (k %% 2 == 0) dif_effect else -dif_effect,
            random = stats::rnorm(1, mean = 0, sd = dif_sd)
          )
        }
        if (dif_par %in% c("a", "both")) {
          drift_a[it, t] <- switch(dif,
            unbalanced = dif_effect_a,
            balanced = if (k %% 2 == 0) dif_effect_a else -dif_effect_a,
            random = stats::rnorm(1, mean = 0, sd = dif_sd_a)
          )
        }
      }
    }
  }

  # shared persons keep their identity and an AR(1) correlated z,
  # z_t = person_cor * z_{t-1} + sqrt(1 - person_cor^2) * e
  n_dep <- round(N * person_overlap)
  ids <- vector("list", n_groups)
  zs <- vector("list", n_groups)
  next_id <- 1L
  for (t in seq_len(n_groups)) {
    if (t == 1L || n_dep == 0L) {
      ids[[t]] <- sprintf("P%06d", seq.int(next_id, next_id + N - 1L))
      next_id <- next_id + N
      # draw z only under dependence, person_overlap = 0 keeps the
      # stream of earlier versions
      if (n_dep > 0L) zs[[t]] <- stats::rnorm(N)
    } else {
      keep <- seq_len(n_dep)
      n_new <- N - n_dep
      new_ids <- if (n_new > 0) {
        sprintf("P%06d", seq.int(next_id, next_id + n_new - 1L))
      } else {
        character(0)
      }
      ids[[t]] <- c(ids[[t - 1L]][keep], new_ids)
      next_id <- next_id + n_new
      z_new <- stats::rnorm(N)
      zs[[t]] <- c(person_cor * zs[[t - 1L]][keep] +
                     sqrt(1 - person_cor^2) * z_new[keep],
                   z_new[-keep])
    }
  }

  data_list <- lapply(seq_len(n_groups), function(t) {
    items <- item_sets[[t]]
    p <- true_ipars[match(items, true_ipars$item), ]
    p$a <- p$a * exp(drift_a[items, t])
    p$b <- p$b + drift[items, t]
    theta <- if (n_dep == 0L) {
      stats::rnorm(N, mean = mu[t], sd = sigma[t])
    } else {
      mu[t] + sigma[t] * zs[[t]]
    }
    if (model == "GPCM") {
      X <- matrix(0L, nrow = N, ncol = length(items))
      for (j in seq_along(items)) {
        thr <- p$b[j] + tau[match(items[j], pool), ]
        P <- em_poly_probs(p$a[j], p$a[j] * cumsum(thr), theta)
        cp <- t(apply(P, 1L, cumsum))
        X[, j] <- as.integer(rowSums(stats::runif(N) > cp))
      }
      storage.mode(X) <- "integer"
    } else {
      prob <- stats::plogis(
        outer(theta, p$b, "-") * matrix(p$a, nrow = N,
                                        ncol = length(items),
                                        byrow = TRUE)
      )
      X <- matrix(stats::rbinom(N * length(items), size = 1, prob = prob),
                  nrow = N)
    }
    colnames(X) <- items
    X <- as.data.frame(X)
    rownames(X) <- ids[[t]]
    X
  })

  structure(data_list,
            true_mu = mu, true_sigma = sigma,
            true_ipars = true_ipars, design = design,
            dif = dif, dif_items = dif_items, dif_effects = drift,
            dif_effects_a = drift_a, person_ids = ids)
}

# item name sets per group for the two designs
build_design_items <- function(n_groups, I, overlap, design) {
  n_common <- round(I * overlap)
  if (n_common < 2) {
    stop("`overlap` must yield at least 2 common items ",
         "(round(I * overlap) >= 2).", call. = FALSE)
  }
  if (design == "all") {
    link_items <- sprintf("A%02d", seq_len(n_common))
    lapply(seq_len(n_groups), function(t) {
      c(link_items, sprintf("W%d_U%02d", t, seq_len(I - n_common)))
    })
  } else {
    # successive, link block L_t shared by groups t and t+1
    if (n_groups > 2 && I - 2 * n_common < 0) {
      stop("the successive design requires I >= 2 * round(I * overlap) so that ",
           "middle groups can hold two link blocks.", call. = FALSE)
    }
    link_blocks <- lapply(seq_len(n_groups - 1), function(t) {
      sprintf("L%d_%02d", t, seq_len(n_common))
    })
    lapply(seq_len(n_groups), function(t) {
      before <- if (t > 1) link_blocks[[t - 1]] else character(0)
      after <- if (t < n_groups) link_blocks[[t]] else character(0)
      n_unique <- I - length(before) - length(after)
      c(before, sprintf("W%d_U%02d", t, seq_len(n_unique)), after)
    })
  }
}
