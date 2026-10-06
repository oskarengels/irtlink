# R/link-weights.R
# quadrature grid and weights for Haebara and Stocking-Lord linking
make_link_weights <- function(weights = c("uniform", "normal_0.5",
                                          "normal_1", "normal_2"),
                              theta = seq(-6, 6, length.out = 101)) {
  weights <- match.arg(weights)
  wgt <- switch(weights,
    uniform = rep(1, length(theta)),
    normal_0.5 = stats::dnorm(theta, mean = 0, sd = 0.5),
    normal_1 = stats::dnorm(theta, mean = 0, sd = 1),
    normal_2 = stats::dnorm(theta, mean = 0, sd = 2)
  )
  list(theta = theta, wgt = wgt)
}
