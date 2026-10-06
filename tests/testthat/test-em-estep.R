test_that("em_estep matches the R reference incl. NA handling", {
  s <- em_ref_setup(seed = 11)
  ref <- em_estep_ref(s$X, s$P, s$pi)
  ed <- em_estep_data(s$X, rep(1, nrow(s$X)))
  out <- em_estep(ed, log(s$P), log(1 - s$P), log(s$pi))
  expect_equal(out$njk, ref$njk, tolerance = 1e-10)
  expect_equal(out$rjk, ref$rjk, tolerance = 1e-10)
  expect_equal(out$loglik, ref$loglik, tolerance = 1e-10)
})

test_that("em_estep loglik equals the direct marginal likelihood", {
  s <- em_ref_setup(seed = 12, N = 20L, J = 4L, K = 11L, na_frac = 0)
  ed <- em_estep_data(s$X, rep(1, nrow(s$X)))
  out <- em_estep(ed, log(s$P), log(1 - s$P), log(s$pi))
  ll <- 0
  for (i in seq_len(nrow(s$X))) {
    lik_k <- s$pi
    for (j in seq_len(ncol(s$X))) {
      lik_k <- lik_k * (if (s$X[i, j] == 1L) s$P[j, ] else 1 - s$P[j, ])
    }
    ll <- ll + log(sum(lik_k))
  }
  expect_equal(out$loglik, ll, tolerance = 1e-10)
})

test_that("weights equal row duplication", {
  s <- em_ref_setup(seed = 13, na_frac = 0)
  X2 <- rbind(s$X, s$X[1:5, ])
  w <- rep(1, nrow(s$X)); w[1:5] <- 2
  dup <- em_estep(em_estep_data(X2, rep(1, nrow(X2))),
                  log(s$P), log(1 - s$P), log(s$pi))
  wtd <- em_estep(em_estep_data(s$X, w),
                  log(s$P), log(1 - s$P), log(s$pi))
  expect_equal(wtd$njk, dup$njk, tolerance = 1e-10)
  expect_equal(wtd$rjk, dup$rjk, tolerance = 1e-10)
  expect_equal(wtd$loglik, dup$loglik, tolerance = 1e-10)
})

test_that("em_estep returns node totals nk", {
  s <- em_ref_setup(seed = 15)
  ed <- em_estep_data(s$X, rep(1, nrow(s$X)))
  out <- em_estep(ed, log(s$P), log(1 - s$P), log(s$pi))
  expect_length(out$nk, length(s$pi))
  expect_equal(sum(out$nk), nrow(s$X), tolerance = 1e-10)
  w <- runif(nrow(s$X), 0.5, 3)
  outw <- em_estep(em_estep_data(s$X, w), log(s$P), log(1 - s$P), log(s$pi))
  expect_equal(sum(outw$nk), sum(w), tolerance = 1e-10)
})

test_that("posterior rows sum to one and njk column sums match", {
  s <- em_ref_setup(seed = 14)
  ed <- em_estep_data(s$X, rep(1, nrow(s$X)))
  out <- em_estep(ed, log(s$P), log(1 - s$P), log(s$pi),
                  return_posterior = TRUE)
  expect_equal(rowSums(out$posterior), rep(1, nrow(s$X)), tolerance = 1e-10)
  expect_equal(rowSums(out$njk), colSums(!is.na(s$X)), tolerance = 1e-10)
})
