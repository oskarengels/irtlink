# tests/testthat/helper-drift-ipars.R
# Identified item parameters with an injected difficulty drift on chosen
# items at groups >= 2 (on the group-1 metric, before re-identification).
# Used to test RMSD detection at the parameter level (no TAM needed).

make_drift_ipars <- function(mu, sigma, I = 20, drift_items = character(0),
                             drift = 0.8) {
  base_a <- rep(c(0.73, 1.25, 1.20, 1.47, 0.97, 1.38, 1.05, 1.14, 1.15, 0.67),
                length.out = I)
  base_b <- rep(c(-1.31, 1.44, -1.20, 0.10, 0.10, -0.74, 1.48, -0.61, 0.82,
                  -0.07), length.out = I)
  item <- sprintf("I%02d", seq_len(I))
  do.call(rbind, lapply(seq_along(mu), function(t) {
    b_t <- base_b
    if (t >= 2) b_t[item %in% drift_items] <- b_t[item %in% drift_items] + drift
    data.frame(
      group = t, item = item,
      a = base_a * sigma[t],
      b = (b_t - mu[t]) / sigma[t]
    )
  }))
}
