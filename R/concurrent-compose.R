# R/concurrent-compose.R
# Concurrent calibration composed across groups for the chain and
# restricted joint approaches. The chain fits successive pairs, the
# restricted joint fits a growing window of groups 1 to t. The joint
# approach itself is a single all-group fit.

# fit a concurrent multiple group model on part of the groups,
# free_items frees those items for partial invariance
concurrent_window <- function(sub_list, model, free_items = NULL,
                              engine = "em", ...) {
  fit_em_concurrent(sub_list, model = model, free_items = free_items, ...)
}

# chain composition, successive-pair fits with accumulated steps
fit_concurrent_chain <- function(data_list, model = "2PL", free_items = NULL,
                                 engine = "em", pweights = NULL, ...) {
  n_groups <- length(data_list)
  trend <- data.frame(group = seq_len(n_groups), mu = NA_real_, sigma = NA_real_)
  trend$mu[1] <- 0; trend$sigma[1] <- 1
  converged <- logical(max(0, n_groups - 1))
  mu_cum <- 0; sigma_cum <- 1
  for (t in seq_len(n_groups - 1)) {
    fit <- concurrent_window(data_list[c(t, t + 1)], model, free_items,
                             engine,
                             pweights = if (is.null(pweights)) NULL
                                        else pweights[c(t, t + 1)], ...)
    mu_step <- fit$trend$mu[2]; sigma_step <- fit$trend$sigma[2]
    mu_cum <- mu_cum + sigma_cum * mu_step
    sigma_cum <- sigma_cum * sigma_step
    trend$mu[t + 1] <- mu_cum; trend$sigma[t + 1] <- sigma_cum
    converged[t] <- isTRUE(fit$converged)
  }
  list(trend = trend, converged = all(converged), ipars = NULL, model = NULL)
}

# restricted joint, growing window over groups 1 to t
fit_concurrent_restricted <- function(data_list, model = "2PL",
                                      free_items = NULL, engine = "em",
                                      pweights = NULL, ...) {
  n_groups <- length(data_list)
  trend <- data.frame(group = seq_len(n_groups), mu = NA_real_, sigma = NA_real_)
  trend$mu[1] <- 0; trend$sigma[1] <- 1
  converged <- logical(max(0, n_groups - 1))
  for (t in seq_len(n_groups - 1L) + 1L) {
    fit <- concurrent_window(data_list[seq_len(t)], model, free_items,
                             engine,
                             pweights = if (is.null(pweights)) NULL
                                        else pweights[seq_len(t)], ...)
    trend$mu[t] <- fit$trend$mu[t]; trend$sigma[t] <- fit$trend$sigma[t]
    converged[t - 1] <- isTRUE(fit$converged)
  }
  list(trend = trend, converged = all(converged), ipars = NULL, model = NULL)
}
