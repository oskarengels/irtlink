# R/rmsd.R

# 2PL item response function on a theta grid
rmsd_irf <- function(theta, a, b) stats::plogis(a * (theta - b))

# parameter-based RMSD per common item per group, pooled reference
# IRF on the group-1 metric, integrated over the group density
rmsd_parameter <- function(ipars, link_method = "mgm", anchor = NULL,
                           theta = seq(-6, 6, length.out = 201)) {
  ipars <- as_ipars(ipars)
  groups <- sort(unique(ipars$group))
  n_groups <- length(groups)
  lk <- link_chain(ipars, method = link_method, anchor = anchor)
  mu <- lk$trend$mu
  sigma <- lk$trend$sigma

  # transformed parameters onto the group-1 reference metric
  ipars$a_star <- ipars$a / sigma[match(ipars$group, groups)]
  ipars$b_star <- mu[match(ipars$group, groups)] +
    sigma[match(ipars$group, groups)] * ipars$b

  # common items present in >= 2 groups
  presence <- table(ipars$item)
  common <- names(presence)[presence >= 2]

  if (length(common) == 0L) {
    return(data.frame(item = character(0), max = numeric(0)))
  }

  # normalized group densities on the grid
  # grid [-6, 6] covers |mu_t| <= 3 and sigma_t <= 1.5
  dens <- lapply(seq_len(n_groups), function(t) {
    w <- stats::dnorm(theta, mean = mu[t], sd = sigma[t])
    w / sum(w)
  })

  # the pooled reference includes the drifted item itself, so use
  # the max column for screening
  out <- lapply(common, function(it) {
    rows <- ipars[ipars$item == it, ]
    wv <- match(rows$group, groups)
    a_bar <- mean(rows$a_star)
    b_bar <- mean(rows$b_star)
    p_ref <- rmsd_irf(theta, a_bar, b_bar)
    rmsd_t <- vapply(seq_along(wv), function(j) {
      p_it <- rmsd_irf(theta, rows$a_star[j], rows$b_star[j])
      sqrt(sum((p_it - p_ref)^2 * dens[[wv[j]]]))
    }, numeric(1))
    res <- stats::setNames(rep(NA_real_, n_groups), paste0("g", groups))
    res[wv] <- rmsd_t
    c(res, max = max(rmsd_t))
  })
  tab <- data.frame(item = common,
                    do.call(rbind, out), row.names = NULL,
                    check.names = FALSE)
  tab
}

# data-based RMSD via a multiple group invariance fit, the _bc
# columns carry the bias-corrected values
rmsd_data <- function(calib, ..., verbose = FALSE) {
  em_rmsd_data(calib)
}
