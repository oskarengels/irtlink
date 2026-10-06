em_multi_sim <- function(seed = 51, n_groups = 3, N = 400L, I = 12L, ...) {
  set.seed(seed)
  sim_trend_data(n_groups = n_groups, N = N, I = I, dif = "none", ...)
}

test_that("em_multi_data builds the union structure", {
  d <- em_multi_sim()
  md <- em_multi_data(d)
  expect_equal(md$n_groups, 3L)
  expect_equal(md$J, length(md$pool))
  expect_equal(md$pool, sort(md$pool))
  for (t in 1:3) {
    wv <- md$groups[[t]]
    expect_equal(md$pool[wv$idx], colnames(as.data.frame(d[[t]])))
    expect_equal(wv$N, nrow(as.data.frame(d[[t]])))
  }
  d2 <- lapply(d, function(x) { x <- as.matrix(x); colnames(x) <- NULL; x })
  expect_error(em_multi_data(d2), "column names")
})

test_that("em_estep_multi reproduces per-group reference results", {
  d <- em_multi_sim(seed = 52)
  md <- em_multi_data(d, collapse = FALSE)
  J <- md$J
  set.seed(1); a <- runif(J, 0.7, 1.6); nu <- rnorm(J)
  g <- matrix(0, J, md$n_groups); g[1, 2] <- 0.4   # one drifted cell
  theta <- em_theta_grid(21L, c(-4, 4))
  pis <- list(dnorm_discrete(theta, 0, 1), dnorm_discrete(theta, 0.3, 1.1),
              dnorm_discrete(theta, 0.6, 1.2))
  es <- em_estep_multi(md, a, nu, g, theta, pis)
  ll_ref <- 0
  for (t in 1:3) {
    X <- as.matrix(as.data.frame(d[[t]])); storage.mode(X) <- "integer"
    idx <- md$groups[[t]]$idx
    P <- em_irf_matrix(a[idx], nu[idx] + g[idx, t], theta)
    ref <- em_estep_ref(X, P, pis[[t]])
    expect_equal(es$groups[[t]]$njk[idx, ], ref$njk, tolerance = 1e-10)
    expect_equal(es$groups[[t]]$rjk[idx, ], ref$rjk, tolerance = 1e-10)
    absent <- setdiff(seq_len(J), idx)
    if (length(absent))
      expect_true(all(es$groups[[t]]$njk[absent, ] == 0))
    expect_equal(sum(es$groups[[t]]$nk), nrow(X), tolerance = 1e-10)
    ll_ref <- ll_ref + ref$loglik
  }
  expect_equal(es$loglik, ll_ref, tolerance = 1e-10)
})

test_that("em_update_normal recovers the moments of the node totals", {
  theta <- em_theta_grid(61L)
  pi_true <- dnorm_discrete(theta, 0.4, 1.2)
  N <- 500
  nk <- N * pi_true
  est <- em_update_normal(nk, theta, N)
  expect_equal(unname(est[1]), sum(theta * pi_true), tolerance = 1e-12)
  expect_equal(unname(est[2]),
               sqrt(sum(theta^2 * pi_true) - sum(theta * pi_true)^2),
               tolerance = 1e-12)
  degen <- numeric(61); degen[31] <- N   # all mass on one node
  est0 <- em_update_normal(degen, theta, N)
  expect_gte(unname(est0[2]), 0.01)
})

em_g_case <- function(seed = 61) {
  set.seed(seed)
  theta <- em_theta_grid(21L, c(-4, 4))
  njk_t <- lapply(1:3, function(t) 100 * dnorm_discrete(theta, 0.2 * t, 1))
  a_true <- 1.2; v_true <- 0.3; g_true <- c(0, 0.5, -0.4)
  rjk_t <- lapply(1:3, function(t) {
    p <- plogis(a_true * theta - v_true - g_true[t])
    njk_t[[t]] * p
  })
  list(theta = theta, njk_t = njk_t, rjk_t = rjk_t,
       a = a_true, nu = v_true, g = g_true)
}

test_that("em_mstep_item_g recovers (a, nu, g) from noiseless counts", {
  cs <- em_g_case()
  res <- em_mstep_item_g(cs$theta, cs$njk_t, cs$rjk_t, tps = 1:3,
                         free_tps = 2:3, a0 = 1, nu0 = 0,
                         g0 = c(0, 0, 0), est_a = TRUE)
  expect_equal(res$a, cs$a, tolerance = 1e-6)
  expect_equal(res$nu, cs$nu, tolerance = 1e-6)
  expect_equal(res$g[2:3], cs$g[2:3], tolerance = 1e-6)
  expect_equal(res$g[1], 0)
})

test_that("em_mstep_item_g matches an optim reference", {
  cs <- em_g_case(62)
  res <- em_mstep_item_g(cs$theta, cs$njk_t, cs$rjk_t, tps = 1:3,
                         free_tps = 2:3, a0 = 1, nu0 = 0,
                         g0 = c(0, 0, 0), est_a = TRUE)
  nll <- function(par) {
    a <- par[1]; nu <- par[2]; g <- c(0, par[3], par[4])
    val <- 0
    for (t in 1:3) {
      p <- pmin(pmax(plogis(a * cs$theta - nu - g[t]), 1e-10), 1 - 1e-10)
      val <- val - sum(cs$rjk_t[[t]] * log(p) +
                       (cs$njk_t[[t]] - cs$rjk_t[[t]]) * log(1 - p))
    }
    val
  }
  o <- stats::optim(c(1, 0, 0, 0), nll, method = "L-BFGS-B",
                    lower = c(0.01, -10, -10, -10),
                    upper = c(10, 10, 10, 10))
  expect_equal(res$a, o$par[1], tolerance = 1e-4)
  expect_equal(res$nu, o$par[2], tolerance = 1e-4)
  expect_equal(res$g[2:3], o$par[3:4], tolerance = 1e-4)
})

test_that("em_mstep_item_g with est_a = FALSE fixes a (1PL)", {
  cs <- em_g_case(63)
  res <- em_mstep_item_g(cs$theta, cs$njk_t, cs$rjk_t, tps = 1:3,
                         free_tps = 2:3, a0 = 1, nu0 = 0,
                         g0 = c(0, 0, 0), est_a = FALSE)
  expect_equal(res$a, 1)
})

test_that("a strong penalty shrinks g to zero, a weak one does not", {
  cs <- em_g_case(65)
  strong <- em_mstep_item_g(cs$theta, cs$njk_t, cs$rjk_t, tps = 1:3,
                            free_tps = 2:3, a0 = 1, nu0 = 0,
                            g0 = c(0, 0, 0), est_a = TRUE,
                            inner_maxit = 200L, inner_tol = 1e-10,
                            pen_logN = 1e6, pen_eps = 0.001)
  expect_lt(max(abs(strong$g[2:3])), 0.05)
  weak <- em_mstep_item_g(cs$theta, cs$njk_t, cs$rjk_t, tps = 1:3,
                          free_tps = 2:3, a0 = 1, nu0 = 0,
                          g0 = c(0, 0, 0), est_a = TRUE,
                          inner_maxit = 200L, inner_tol = 1e-10,
                          pen_logN = 0.01, pen_eps = 0.001)
  expect_equal(weak$g[2:3], cs$g[2:3], tolerance = 0.05)
})

test_that("the penalized update matches an optim reference", {
  cs <- em_g_case(66)
  logN <- log(5000); eps <- 0.01
  res <- em_mstep_item_g(cs$theta, cs$njk_t, cs$rjk_t, tps = 1:3,
                         free_tps = 2:3, a0 = 1, nu0 = 0,
                         g0 = c(0, 0, 0), est_a = TRUE,
                         inner_maxit = 500L, inner_tol = 1e-12,
                         pen_logN = logN, pen_eps = eps)
  nqf <- function(par) {
    a <- par[1]; nu <- par[2]; g <- c(0, par[3], par[4])
    val <- 0
    for (t in 1:3) {
      p <- pmin(pmax(plogis(a * cs$theta - nu - g[t]), 1e-10), 1 - 1e-10)
      val <- val - sum(cs$rjk_t[[t]] * log(p) +
                       (cs$njk_t[[t]] - cs$rjk_t[[t]]) * log(1 - p))
    }
    val + logN * sum(g[2:3]^2 / (g[2:3]^2 + eps))
  }
  o <- stats::optim(c(res$a, res$nu, res$g[2], res$g[3]), nqf,
                    method = "L-BFGS-B",
                    lower = c(0.01, -10, -10, -10),
                    upper = c(10, 10, 10, 10))
  # same or better penalized objective than optim started at our solution
  ours <- nqf(c(res$a, res$nu, res$g[2], res$g[3]))
  expect_lt(ours - o$value, 1e-3)
})

test_that("the unpenalized path is unchanged", {
  cs <- em_g_case(62)
  r0 <- em_mstep_item_g(cs$theta, cs$njk_t, cs$rjk_t, tps = 1:3,
                        free_tps = 2:3, a0 = 1, nu0 = 0,
                        g0 = c(0, 0, 0), est_a = TRUE)
  r1 <- em_mstep_item_g(cs$theta, cs$njk_t, cs$rjk_t, tps = 1:3,
                        free_tps = 2:3, a0 = 1, nu0 = 0,
                        g0 = c(0, 0, 0), est_a = TRUE,
                        pen_logN = 0, pen_eps = NULL)
  expect_identical(r0, r1)
})

test_that("em_parmap_concurrent frees h for slope-freed items", {
  d <- em_multi_sim(seed = 53)
  md <- em_multi_data(d)
  pm <- em_parmap_concurrent(md, "2PL", free_items = NULL,
                             free_slope_items = md$pool[2])
  hh <- pm[pm$parname == "h", ]
  expect_true(all(hh$est[hh$item == md$pool[2] & hh$group >= 2L]))
  expect_true(all(!hh$est[hh$item != md$pool[2]]))
  expect_true(all(!hh$est[hh$group == 1L]))
})

em_h_case <- function(seed = 67) {
  set.seed(seed)
  theta <- em_theta_grid(21L, c(-4, 4))
  njk_t <- lapply(1:3, function(t) 100 * dnorm_discrete(theta, 0.2 * t, 1))
  a_true <- 1.1; v_true <- 0.2
  g_true <- c(0, 0.4, 0); h_true <- c(0, -0.3, 0.35)
  rjk_t <- lapply(1:3, function(t) {
    p <- plogis((a_true + h_true[t]) * theta - v_true - g_true[t])
    njk_t[[t]] * p
  })
  list(theta = theta, njk_t = njk_t, rjk_t = rjk_t,
       a = a_true, nu = v_true, g = g_true, h = h_true)
}

test_that("em_mstep_item_g recovers slope offsets from noiseless counts", {
  cs <- em_h_case()
  res <- em_mstep_item_g(cs$theta, cs$njk_t, cs$rjk_t, tps = 1:3,
                         free_tps = 2:3, a0 = 1, nu0 = 0,
                         g0 = c(0, 0, 0), est_a = TRUE,
                         inner_maxit = 200L, inner_tol = 1e-12,
                         h0 = c(0, 0, 0), free_h_tps = 2:3)
  expect_equal(res$a, cs$a, tolerance = 1e-5)
  expect_equal(res$nu, cs$nu, tolerance = 1e-5)
  expect_equal(res$g[2:3], cs$g[2:3], tolerance = 1e-5)
  expect_equal(res$h[2:3], cs$h[2:3], tolerance = 1e-5)
  expect_equal(res$h[1], 0)
})

test_that("slope-only freeing recovers h with g fixed at zero", {
  cs <- em_h_case(68)
  # data without intercept drift so g = 0 is the truth
  rjk_t <- lapply(1:3, function(t) {
    p <- plogis((cs$a + cs$h[t]) * cs$theta - cs$nu)
    cs$njk_t[[t]] * p
  })
  res <- em_mstep_item_g(cs$theta, cs$njk_t, rjk_t, tps = 1:3,
                         free_tps = integer(0), a0 = 1, nu0 = 0,
                         g0 = c(0, 0, 0), est_a = TRUE,
                         inner_maxit = 200L, inner_tol = 1e-12,
                         h0 = c(0, 0, 0), free_h_tps = 2:3)
  expect_equal(res$h[2:3], cs$h[2:3], tolerance = 1e-5)
  expect_equal(res$g, c(0, 0, 0))
})

test_that("h defaults keep the Stage-3 behavior identical", {
  cs <- em_g_case(62)
  r0 <- em_mstep_item_g(cs$theta, cs$njk_t, cs$rjk_t, tps = 1:3,
                        free_tps = 2:3, a0 = 1, nu0 = 0,
                        g0 = c(0, 0, 0), est_a = TRUE)
  expect_equal(r0$h, c(0, 0, 0))
  expect_equal(r0$g[2:3], cs$g[2:3], tolerance = 1e-4)
})

test_that("pooled kernel equals tied multiple group update for invariant items", {
  cs <- em_g_case(64)
  npool <- Reduce(`+`, cs$njk_t); rpool <- Reduce(`+`, cs$rjk_t)
  cpp <- em_mstep_items_cpp(cs$theta, matrix(npool, 1), matrix(rpool, 1),
                            1, 0, TRUE, 0.01, 10, -10, 10, 100L, 1e-10)
  rg <- em_mstep_item_g(cs$theta, cs$njk_t, cs$rjk_t, tps = 1:3,
                        free_tps = integer(0), a0 = 1, nu0 = 0,
                        g0 = c(0, 0, 0), est_a = TRUE,
                        inner_maxit = 100L, inner_tol = 1e-12)
  expect_equal(rg$a, cpp$a[1], tolerance = 1e-6)
  expect_equal(rg$nu, cpp$nu[1], tolerance = 1e-6)
})

test_that("em_parmap_concurrent ties a/nu and frees g correctly", {
  d <- em_multi_sim()
  md <- em_multi_data(d)
  it_free <- md$pool[1]
  pm <- em_parmap_concurrent(md, "2PL", free_items = it_free)
  av <- pm[pm$parname %in% c("a", "nu"), ]
  for (it in md$pool) {
    for (pn in c("a", "nu")) {
      px <- unique(av$parindex[av$item == it & av$parname == pn])
      expect_length(px, 1L)   # one shared parindex per item per parameter
    }
  }
  gg <- pm[pm$parname == "g", ]
  expect_true(all(!gg$est[gg$group == 1L]))
  expect_true(all(gg$est[gg$item == it_free & gg$group >= 2L]))
  expect_true(all(!gg$est[gg$item != it_free]))
  expect_true(all(gg$value[!gg$est] == 0))
  pm1 <- em_parmap_concurrent(md, "1PL", free_items = NULL)
  expect_true(all(!pm1$est[pm1$parname == "a"]))
  expect_true(all(!pm1$est[pm1$parname == "g"]))
  expect_warning(em_parmap_concurrent(md, "2PL", free_items = "nope"),
                 "not found")
})
