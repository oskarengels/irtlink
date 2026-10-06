# generalized partial credit model calibration and linking

ip_gpcm <- function(I = 12, K = 2) {
  b_ref <- seq(-1.2, 1.2, length.out = I)
  a_ref <- rep(c(0.8, 1.2, 1.0), length.out = I)
  tau_v <- seq(-0.5, 0.5, length.out = K)
  ip <- data.frame(
    group = rep(1:2, each = I),
    item = rep(sprintf("I%02d", seq_len(I)), 2),
    a = c(a_ref, a_ref * 1.1),
    b = c(b_ref, (b_ref - 0.3) / 1.1)
  )
  for (v in seq_len(K)) {
    ip[[paste0("tau", v)]] <- c(rep(tau_v[v], I), rep(tau_v[v] / 1.1, I))
  }
  ip
}

sim_gpcm <- function(n_groups = 2, N = 800, seed = 7) {
  set.seed(seed)
  sim_trend_data(n_groups = n_groups, N = N, I = 10, overlap = 1,
                 mu = 0.3 * (seq_len(n_groups) - 1),
                 sigma = rep(1, n_groups), model = "GPCM", n_cat = 3)
}

test_that("every chain method recovers an exact GPCM transformation", {
  ip <- ip_gpcm()
  for (m in c("mgm", "mm", "haberman", "haebara", "sl")) {
    lk <- link_groups(ip, approach = "chain", method = m)
    expect_equal(lk$trend$mu[2], 0.3, tolerance = 1e-5)
    expect_equal(lk$trend$sigma[2], 1.1, tolerance = 1e-5)
  }
  lk <- link_groups(ip, approach = "joint", method = "haberman")
  expect_equal(lk$trend$mu[2], 0.3, tolerance = 1e-5)
})

test_that("GPCM calibration recovers simulated parameters", {
  d <- sim_gpcm()
  cal <- calibrate(d, model = "GPCM")
  expect_true(all(cal$converged))
  expect_true(all(c("tau1", "tau2") %in% names(cal$ipars)))
  tru <- attr(d, "true_ipars")
  g1 <- cal$ipars[cal$ipars$group == 1, ]
  m <- match(g1$item, tru$item)
  expect_gt(cor(g1$a, tru$a[m]), 0.7)
  expect_lt(mean(abs(g1$b - tru$b[m])), 0.15)
})

test_that("PCM fixes all discriminations at one", {
  d <- sim_gpcm(N = 300)
  cal <- calibrate(d, model = "PCM")
  expect_true(all(cal$ipars$a == 1))
})

test_that("GPCM on dichotomous data matches the 2PL fit", {
  set.seed(11)
  d <- sim_trend_data(n_groups = 2, N = 1000, I = 10, overlap = 1)
  c2 <- calibrate(d, model = "2PL")
  cg <- calibrate(d, model = "GPCM")
  expect_lt(max(abs(c2$ipars$a - cg$ipars$a)), 0.02)
  expect_lt(max(abs(c2$ipars$b - cg$ipars$b)), 0.02)
})

test_that("jackknife and ajk agree for a GPCM chain link", {
  ip <- ip_gpcm()
  ip$b[ip$group == 2][2] <- ip$b[ip$group == 2][2] + 0.5
  lk <- link_groups(ip, approach = "chain", method = "mgm")
  jk <- linking_error(lk, method = "jackknife")
  aj <- linking_error(lk, method = "ajk")
  expect_gt(jk$le$le_mu[2], 0)
  expect_equal(aj$le$le_mu[2], jk$le$le_mu[2], tolerance = 0.1)
})

test_that("pairwise joint criteria recover an exact GPCM setting", {
  I <- 10
  K <- 2
  mu_t <- c(0, 0.3, 0.5)
  sigma_t <- c(1, 1.1, 1.2)
  b_ref <- seq(-1.2, 1.2, length.out = I)
  a_ref <- rep(c(0.8, 1.2, 1.0), length.out = I)
  tau_v <- c(-0.5, 0.5)
  ip <- do.call(rbind, lapply(1:3, function(t) {
    d <- data.frame(group = t, item = sprintf("I%02d", seq_len(I)),
                    a = a_ref * sigma_t[t],
                    b = (b_ref - mu_t[t]) / sigma_t[t])
    for (v in seq_len(K)) d[[paste0("tau", v)]] <- tau_v[v] / sigma_t[t]
    d
  }))
  for (m in c("haebara", "sl")) {
    lk <- link_groups(ip, approach = "joint", method = m,
                      variant = "pairwise")
    expect_equal(lk$trend$mu, mu_t, tolerance = 1e-3)
    expect_equal(lk$trend$sigma, sigma_t, tolerance = 1e-3)
  }
})

test_that("polytomous guards fire", {
  ip <- ip_gpcm()
  expect_error(link_groups(ip, approach = "joint", method = "haebara"),
               "polytomous")
  lk <- link_groups(ip, approach = "chain", method = "mgm")
  expect_error(linking_error(lk, method = "sandwich_esw"), "polytomous")
  lh <- link_groups(ip, approach = "chain", method = "haebara")
  expect_error(linking_error(lh, method = "ajk"), "polytomous")
  d <- sim_gpcm(N = 200)
  expect_error(calibrate(d, model = "GPCM", calibration = "concurrent"),
               "separate")
  cal <- calibrate(d, model = "GPCM")
  expect_error(detect_dif(cal), "polytomous")
  expect_error(person_scores(cal, method = "wle"), "eap")
})

test_that("category gaps are rejected", {
  d <- sim_gpcm(N = 200)
  d[[1]][[1]][d[[1]][[1]] == 1] <- 2
  expect_error(calibrate(d, model = "GPCM"), "skips response categories")
})

test_that("the GPCM simulator returns categories and attributes", {
  d <- sim_gpcm(N = 150)
  vals <- sort(unique(unlist(d)))
  expect_identical(vals, c(0L, 1L, 2L))
  expect_true("tau1" %in% names(attr(d, "true_ipars")))
  expect_error(sim_trend_data(n_groups = 2, N = 50, model = "GPCM",
                              n_cat = 1), "n_cat")
})

test_that("estimate_trend runs the GPCM separate pipeline", {
  d <- sim_gpcm(N = 500)
  tr <- estimate_trend(d, link = "mgm", model = "GPCM")
  expect_equal(tr$trend$mu[2], 0.3, tolerance = 0.15)
  expect_error(estimate_trend(d, model = "GPCM", dif = "purify"),
               "polytomous")
})
