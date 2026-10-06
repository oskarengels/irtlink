link_chain <- function(x, method = c("mgm", "mm", "haberman", "haebara", "sl"),
                       pow = 2, weights = "uniform", use_intercepts = TRUE,
                       anchor = NULL, type = c("asymm", "symm"),
                       theta_metric = c("to", "from"), ref = NULL, ...) {
  method <- match.arg(method)
  type <- match.arg(type)
  theta_metric <- match.arg(theta_metric)
  if ((type != "asymm" || theta_metric != "to") &&
      !method %in% c("haebara", "sl")) {
    stop("`type` and `theta_metric` apply to method = \"haebara\" or ",
         "\"sl\" only.", call. = FALSE)
  }
  ipars <- as_ipars(x)
  groups <- sort(unique(ipars$group))
  n_groups <- length(groups)
  if (n_groups < 2) {
    stop("Chain linking requires at least two groups.", call. = FALSE)
  }
  estimator <- switch(method,
    mm = est_link_mm,
    mgm = est_link_mgm,
    haberman = est_link_haberman,
    haebara = est_link_haebara,
    sl = est_link_sl
  )
  est_args <- switch(method,
    haberman = list(pow = pow, use_intercepts = use_intercepts),
    haebara = ,
    sl = list(pow = pow, weights = weights, type = type,
              theta_metric = theta_metric),
    list()  # mm/mgm: closed form; pow/weights/use_intercepts are ignored
  )

  trend <- data.frame(group = groups, mu = NA_real_, sigma = NA_real_)
  trend$mu[1] <- 0
  trend$sigma[1] <- 1
  steps <- vector("list", n_groups - 1)
  mu_cum <- 0
  sigma_cum <- 1
  for (t in seq_len(n_groups - 1)) {
    pars <- ipars_pair(ipars, groups[t], groups[t + 1])
    if (!is.null(anchor)) {
      pars <- pars[pars$item %in% anchor, , drop = FALSE]
      if (nrow(pars) == 0) {
        stop(sprintf(
          "No anchor items among the common items of groups %d and %d.",
          groups[t], groups[t + 1]), call. = FALSE)
      }
      if (nrow(pars) < length(anchor)) {
        message(sprintf(
          "Group %d->%d: %d of %d anchor item(s) are common to this pair and will be used.",
          groups[t], groups[t + 1], nrow(pars), length(anchor)))
      }
    }
    est <- do.call(estimator, c(list(pars), est_args, list(...)))
    mu_cum <- mu_cum + sigma_cum * est$mu
    sigma_cum <- sigma_cum * est$sigma
    trend$mu[t + 1] <- mu_cum
    trend$sigma[t + 1] <- sigma_cum
    steps[[t]] <- data.frame(
      from = groups[t], to = groups[t + 1],
      mu_step = est$mu, sigma_step = est$sigma,
      n_common = nrow(pars), converged = est$converged
    )
  }

  trend <- rebase_trend(trend, ref)

  structure(
    list(
      method = method,
      approach = "chain",
      trend = trend,
      steps = do.call(rbind, steps),
      ipars = ipars,
      anchor = anchor,
      ref = ref,
      control = c(list(pow = pow, weights = weights,
                       use_intercepts = use_intercepts, type = type,
                       theta_metric = theta_metric, ref = ref), list(...)),
      le = NULL,
      call = match.call()
    ),
    class = "irtlink_link"
  )
}

# re-express the trend in the metric of the reference group,
# mu' = (mu - mu_ref) / sigma_ref and sigma' = sigma / sigma_ref
rebase_trend <- function(trend, ref) {
  if (is.null(ref)) return(trend)
  r <- match(ref, trend$group)
  if (is.na(r)) {
    stop("`ref` must be one of the linked groups (",
         paste(trend$group, collapse = ", "), "); got ", ref, ".",
         call. = FALSE)
  }
  if (r == 1L) return(trend)
  mu_r <- trend$mu[r]
  sigma_r <- trend$sigma[r]
  trend$mu <- (trend$mu - mu_r) / sigma_r
  trend$sigma <- trend$sigma / sigma_r
  trend
}
