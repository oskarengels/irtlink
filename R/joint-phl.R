# R/joint-phl.R
# pairwise Haberman linking (PHL) over all group pairs
# (Robitzsch, 2025, AppliedMath), two closed-form weighted least
# squares stages
joint_phl <- function(ipars, item_weights = c("inverse_admin", "uniform")) {
  variant <- match.arg(item_weights)
  groups <- sort(unique(ipars$group))
  Tn <- length(groups)
  if (Tn < 2) {
    stop("PHL linking requires at least two groups.", call. = FALSE)
  }
  by_group <- split(ipars, ipars$group)
  items_t <- lapply(by_group, function(d) d$item)
  aL <- lapply(by_group, function(d) stats::setNames(d$a, d$item))
  bL <- lapply(by_group, function(d) stats::setNames(d$b, d$item))

  all_a <- unlist(aL, use.names = FALSE)
  if (any(all_a <= 0)) {
    stop("Non-positive discrimination(s) in item parameters; ",
         "remove degenerate items before PHL linking.", call. = FALSE)
  }

  presence <- table(ipars$item)  # G_i per item
  omega <- if (variant == "uniform") {
    stats::setNames(rep(1, length(presence)), names(presence))
  } else {
    stats::setNames(as.numeric(1 / presence), names(presence))
  }

  P <- Tn - 1
  idx <- function(t) t - 1  # group t (>= 2) -> parameter column

  # stage 1, scale parameters s with sigma = exp(s)
  XtWX_s <- matrix(0, P, P)
  XtWy_s <- numeric(P)
  rows <- list()
  r <- 0L
  for (g in 1:(Tn - 1)) {
    for (h in (g + 1):Tn) {
      com <- intersect(items_t[[g]], items_t[[h]])
      for (it in com) {
        w <- omega[[it]]
        ys <- log(aL[[g]][[it]]) - log(aL[[h]][[it]])
        x <- numeric(P)
        if (g >= 2) x[idx(g)] <- x[idx(g)] + 1
        if (h >= 2) x[idx(h)] <- x[idx(h)] - 1
        XtWX_s <- XtWX_s + w * (x %o% x)
        XtWy_s <- XtWy_s + w * x * ys
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
           "disconnected (no item path between some pairs of groups). ",
           "Original error: ", conditionMessage(e), call. = FALSE)
    }
  )
  s <- c(0, s_free)
  sigma <- exp(s)

  # stage 2, means on the stage-1 scales
  XtWX_m <- XtWX_s
  XtWy_m <- numeric(P)
  for (rr in rows) {
    ym <- -(sigma[rr$g] * bL[[rr$g]][[rr$it]] -
            sigma[rr$h] * bL[[rr$h]][[rr$it]])
    XtWy_m <- XtWy_m + rr$w * rr$x * ym
  }
  mu_free <- tryCatch(
    solve(XtWX_m, XtWy_m),
    error = function(e) {
      stop("PHL mean system is singular; the group link graph may be ",
           "disconnected (no item path between some pairs of groups). ",
           "Original error: ", conditionMessage(e), call. = FALSE)
    }
  )
  mu <- c(0, mu_free)

  list(mu = unname(mu), sigma = unname(sigma), converged = TRUE)
}
