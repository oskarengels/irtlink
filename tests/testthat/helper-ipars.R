# tests/testthat/helper-ipars.R
# Identified (infinite-N) item parameters for n_groups with known mu/sigma
# and no drift: a_hat = a * sigma_t, b_hat = (b - mu_t) / sigma_t.
# Linking recovers mu/sigma exactly, so tests can use tight tolerances.

make_identified_ipars <- function(mu, sigma, I = 20) {
  base_a <- rep(c(0.73, 1.25, 1.20, 1.47, 0.97, 1.38, 1.05, 1.14, 1.15, 0.67),
                length.out = I)
  base_b <- rep(c(-1.31, 1.44, -1.20, 0.10, 0.10, -0.74, 1.48, -0.61, 0.82,
                  -0.07), length.out = I)
  do.call(rbind, lapply(seq_along(mu), function(t) {
    data.frame(
      group = t,
      item = sprintf("I%02d", seq_len(I)),
      a = base_a * sigma[t],
      b = (base_b - mu[t]) / sigma[t]
    )
  }))
}
