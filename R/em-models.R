# R/em-models.R
# item model registry for the internal EM engine, each model defines
# parameter names, start values, estimation flags, and bounds

em_model_registry <- function(model = c("2PL", "1PL", "3PL")) {
  model <- match.arg(model)
  if (model == "3PL") {
    return(list(
      id = 3L,
      parnames = c("a", "nu", "c"),
      start = c(a = 1, nu = 0, c = 0.15),
      est = c(a = TRUE, nu = TRUE, c = TRUE),
      lower = c(a = 0.1, nu = -10, c = 0.001),
      upper = c(a = 10, nu = 10, c = 0.5)
    ))
  }
  list(
    id = if (model == "2PL") 2L else 1L,
    parnames = c("a", "nu"),
    start = c(a = 1, nu = 0),
    est = c(a = model == "2PL", nu = TRUE),
    lower = c(a = 0.1, nu = -10),
    upper = c(a = 10, nu = 10)
  )
}

# 3PL item M step, Bayes modal with a beta prior on the guessing
# parameter (Hanson, 1998, Eqs. 29-30)
em_mstep_item_3pl_objective <- function(par, theta, njk, rjk, prior_c) {
  a <- par[1]; nu <- par[2]; cc <- par[3]
  psi <- stats::plogis(a * theta - nu)
  P <- pmin(pmax(cc + (1 - cc) * psi, 1e-10), 1 - 1e-10)
  sum(rjk * log(P) + (njk - rjk) * log(1 - P)) +
    (prior_c[1] - 1) * log(cc) + (prior_c[2] - 1) * log(1 - cc)
}

em_mstep_item_3pl_gradient <- function(par, theta, njk, rjk, prior_c) {
  a <- par[1]; nu <- par[2]; cc <- par[3]
  psi <- stats::plogis(a * theta - nu)
  P <- pmin(pmax(cc + (1 - cc) * psi, 1e-10), 1 - 1e-10)
  w <- rjk / P - (njk - rjk) / (1 - P)
  dPdeta <- (1 - cc) * psi * (1 - psi)
  c(sum(w * dPdeta * theta),
    -sum(w * dPdeta),
    sum(w * (1 - psi)) +
      (prior_c[1] - 1) / cc - (prior_c[2] - 1) / (1 - cc))
}

em_mstep_item_3pl <- function(theta, njk, rjk, a0, nu0, c0,
                              prior_c = c(5, 17),
                              lower_a = 0.1, upper_a = 10) {
  o <- stats::optim(c(a0, nu0, c0),
                    fn = function(p) -em_mstep_item_3pl_objective(
                      p, theta, njk, rjk, prior_c),
                    gr = function(p) -em_mstep_item_3pl_gradient(
                      p, theta, njk, rjk, prior_c),
                    method = "L-BFGS-B",
                    lower = c(lower_a, -10, 0.001),
                    upper = c(upper_a, 10, 0.5))
  list(a = o$par[1], nu = o$par[2], c = o$par[3])
}
