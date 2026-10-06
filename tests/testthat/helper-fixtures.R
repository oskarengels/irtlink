# shared test fixtures

make_drift_calib <- function(drift_items, drift = 1.0, mu = c(0, 0.3, 0.6),
                             sigma = c(1, 1.1, 1.2), I = 20) {
  ip <- as_ipars(make_drift_ipars(mu, sigma, I = I,
                                  drift_items = drift_items, drift = drift))
  structure(list(ipars = ip, model = "2PL", calibration = "separate",
                 engine = "em", n_groups = length(mu), models = NULL),
            class = "irtlink_calib")
}

# identified parameters with random drift at the later groups
drifted_ipars <- function(n_groups, seed, drift = 1,
                          sd_b = 0.25, sd_log_a = 0.15) {
  mu <- c(0, 0.35, 0.6)[seq_len(n_groups)]
  sigma <- c(1, 1.15, 1.25)[seq_len(n_groups)]
  ip <- make_identified_ipars(mu, sigma, I = 12)
  set.seed(seed)
  for (t in seq_len(n_groups)[-1]) {
    sel <- ip$group == t
    ip$b[sel] <- ip$b[sel] + stats::rnorm(sum(sel), sd = drift * sd_b)
    ip$a[sel] <- ip$a[sel] * exp(stats::rnorm(sum(sel),
                                              sd = drift * sd_log_a))
  }
  ip
}

ajk_hae_ipars <- function(n_groups = 2, seed = 11, drift = 1) {
  drifted_ipars(n_groups, seed, drift, sd_b = 0.25, sd_log_a = 0.15)
}

jh_ipars <- function(n_groups = 3, seed = 21, drift = 1) {
  drifted_ipars(n_groups, seed, drift, sd_b = 0.2, sd_log_a = 0.1)
}

ajk_sl_ipars <- ajk_hae_ipars
jp_ipars <- jh_ipars
