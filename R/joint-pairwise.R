# R/joint-pairwise.R
# pairwise Haebara and Stocking-Lord linking over all group pairs,
# item response functions (Haebara) and test response functions
# (Stocking-Lord), for polytomous items the category response and
# expected score functions (Kim & Lee, 2006)

est_joint_pairwise <- function(ipars, method = c("Hae", "SL"), theta,
                               weights = "normal_1", maxit = 1000L) {
  method <- match.arg(method)
  w <- make_link_weights(weights, theta)
  qw <- w$wgt / sum(w$wgt)
  th <- w$theta
  poly <- ipars_poly(ipars)
  tau_cols <- grep("^tau[0-9]+$", names(ipars), value = TRUE)
  groups <- sort(unique(ipars$group))
  NM <- length(groups)
  by_group <- lapply(groups, function(wv)
    ipars[ipars$group == wv, , drop = FALSE])
  # thresholds b + tau_v of the common items of one group
  thr_of <- function(d, com) {
    tau <- as.matrix(d[match(com, d$item), tau_cols, drop = FALSE])
    lapply(seq_along(com), function(i) {
      tv <- tau[i, !is.na(tau[i, ])]
      if (!length(tv)) tv <- 0
      d$b[match(com[i], d$item)] + tv
    })
  }
  combis <- utils::combn(NM, 2L)
  pairs <- list()
  for (j in seq_len(ncol(combis))) {
    p <- combis[1L, j]
    q <- combis[2L, j]
    com <- intersect(by_group[[p]]$item, by_group[[q]]$item)
    if (!length(com)) next
    dp <- by_group[[p]]
    dq <- by_group[[q]]
    cp <- if (is.null(dp$c)) rep(0, length(com)) else dp$c[match(com, dp$item)]
    cq <- if (is.null(dq$c)) rep(0, length(com)) else dq$c[match(com, dq$item)]
    pr <- list(
      p = p, q = q,
      ap = dp$a[match(com, dp$item)], bp = dp$b[match(com, dp$item)],
      aq = dq$a[match(com, dq$item)], bq = dq$b[match(com, dq$item)],
      cp = cp, cq = cq
    )
    if (poly) {
      pr$thrp <- thr_of(dp, com)
      pr$thrq <- thr_of(dq, com)
      Kp <- lengths(pr$thrp)
      if (!identical(Kp, lengths(pr$thrq))) {
        stop("Common item(s) between groups ", groups[p], " and ",
             groups[q], " have a different number of observed ",
             "categories.", call. = FALSE)
      }
      pr$K <- Kp
    }
    pairs[[length(pairs) + 1L]] <- pr
  }
  if (!length(pairs)) {
    stop("Pairwise joint linking found no pair of groups with common ",
         "items.", call. = FALSE)
  }
  cat_probs <- function(a_j, thr_j, tg) {
    em_poly_probs(a_j, a_j * cumsum(thr_j), tg)
  }
  crit <- function(x) {
    bvec <- c(0, x[seq_len(NM - 1L)])
    avec <- c(1, x[NM - 1L + seq_len(NM - 1L)])
    val <- 0
    for (pr in pairs) {
      th_p <- avec[pr$p] * th + bvec[pr$p]
      th_q <- avec[pr$q] * th + bvec[pr$q]
      if (poly) {
        if (method == "Hae") {
          for (i in seq_along(pr$ap)) {
            e <- cat_probs(pr$ap[i], pr$thrp[[i]], th_p)[, -1,
                                                         drop = FALSE] -
              cat_probs(pr$aq[i], pr$thrq[[i]], th_q)[, -1, drop = FALSE]
            val <- val + sum(qw * rowSums(e^2))
          }
        } else {
          T1 <- 0
          T2 <- 0
          for (i in seq_along(pr$ap)) {
            T1 <- T1 + as.numeric(cat_probs(pr$ap[i], pr$thrp[[i]],
                                            th_p) %*% (0:pr$K[i]))
            T2 <- T2 + as.numeric(cat_probs(pr$aq[i], pr$thrq[[i]],
                                            th_q) %*% (0:pr$K[i]))
          }
          val <- val + sum(qw * (T1 - T2)^2)
        }
        next
      }
      P1 <- pr$cp + (1 - pr$cp) *
        stats::plogis(outer(pr$ap, th_p) - pr$ap * pr$bp)
      P2 <- pr$cq + (1 - pr$cq) *
        stats::plogis(outer(pr$aq, th_q) - pr$aq * pr$bq)
      if (method == "Hae") {
        val <- val + sum((P1 - P2)^2 %*% qw)
      } else {
        d <- colSums(P1) - colSums(P2)
        val <- val + sum(qw * d^2)
      }
    }
    val
  }
  opt <- stats::optim(c(rep(0, NM - 1L), rep(1, NM - 1L)), crit,
                      method = "L-BFGS-B",
                      control = list(maxit = maxit, factr = 1e4))
  bvec <- c(0, opt$par[seq_len(NM - 1L)])
  avec <- c(1, opt$par[NM - 1L + seq_len(NM - 1L)])
  list(mu = -bvec / avec, sigma = 1 / avec,
       converged = opt$convergence == 0, joint_fit = NULL)
}
