test_that("person scores sum to the analytic marginal score", {
  set.seed(431)
  d <- sim_trend_data(n_groups = 2, N = 400, I = 8, dif = "none")
  fit <- em_fit(as.matrix(as.data.frame(d[[1]])), "2PL")
  S <- em_person_scores(fit)
  expect_equal(nrow(S), nrow(fit$dat))
  expect_equal(ncol(S), 2L * length(fit$item))
  expect_equal(as.numeric(colSums(S)),
               em_score_ipars(fit, fit$a, fit$nu), tolerance = 1e-8)
  # at the EM solution the total gradient is near zero
  expect_lt(max(abs(colSums(S))) / nrow(S), 1e-3)
})

test_that("person scores handle missing responses (booklets)", {
  set.seed(432)
  d <- sim_trend_data(n_groups = 2, N = 400, I = 8, dif = "none")
  X <- as.matrix(as.data.frame(d[[1]]))
  X[sample(length(X), floor(0.3 * length(X)))] <- NA
  fit <- em_fit(X, "2PL")
  S <- em_person_scores(fit)
  # a person contributes zero score to items not administered
  miss <- which(is.na(X[1, ]))
  expect_true(length(miss) == 0 ||
                all(S[1, as.vector(rbind(2 * miss - 1, 2 * miss))] == 0))
  expect_equal(as.numeric(colSums(S)),
               em_score_ipars(fit, fit$a, fit$nu), tolerance = 1e-8)
})

test_that("stacked_vcov: structure, independence block, and cross blocks", {
  set.seed(433)
  d <- sim_trend_data(n_groups = 2, N = 500, I = 8, dif = "none")
  d <- lapply(d, function(X) as.matrix(as.data.frame(X)))
  cal <- suppressWarnings(calibrate(d, model = "2PL", calibration = "separate",
                   engine = "em", keep_vcov = TRUE))
  # disjoint persons: no cross-group meat, vcov equals vcov_nocross
  ids_disjoint <- list(paste0("p", 1:500), paste0("q", 1:500))
  sv0 <- stacked_vcov(cal, ids_disjoint)
  expect_s3_class(sv0, "irtlink_stacked_vcov")
  expect_equal(sv0$vcov, sv0$vcov_nocross)
  expect_equal(unname(sv0$overlap), matrix(c(500L, 0L, 0L, 500L), 2))
  # the independence matrix reproduces the per-group vcov_ipars exactly
  for (t in 1:2) {
    nm <- rownames(sv0$vcov)[startsWith(rownames(sv0$vcov),
                                        paste0("g", t, ":"))]
    Vt <- sv0$vcov_independence[nm, nm]
    dimnames(Vt) <- dimnames(cal$vcov_ipars[[t]])
    expect_equal(Vt, cal$vcov_ipars[[t]], tolerance = 1e-10)
  }
  # cross-group independence block is exactly zero
  w1 <- startsWith(rownames(sv0$vcov), "g1:")
  expect_true(all(sv0$vcov_independence[w1, !w1] == 0))
  expect_true(all(sv0$vcov[w1, !w1] == 0))
  # overlapping persons produce nonzero cross blocks and keep symmetry
  ids_overlap <- list(paste0("p", 1:500), paste0("p", 201:700))
  sv1 <- stacked_vcov(cal, ids_overlap)
  expect_equal(sv1$overlap["g1", "g2"], 300L)
  expect_gt(max(abs(sv1$vcov[w1, !w1])), 0)
  expect_equal(sv1$vcov, t(sv1$vcov))
  # per-item blocks carry both groups of a common item
  it <- cal$ipars$item[1]
  expect_equal(nrow(sv1$vcov_items[[it]]), 4L)
})

test_that("stacked_vcov feeds linking_error like the per-group list", {
  set.seed(434)
  d <- sim_trend_data(n_groups = 2, N = 500, I = 8, dif = "none")
  d <- lapply(d, function(X) as.matrix(as.data.frame(X)))
  cal <- suppressWarnings(calibrate(d, model = "2PL", calibration = "separate",
                   engine = "em", keep_vcov = TRUE))
  link <- link_chain(cal, method = "mgm")
  le_list <- linking_error(link, method = "jackknife_bc",
                           vcov = cal$vcov_ipars)
  ids <- list(paste0("p", 1:500), paste0("p", 1:500))
  sv <- stacked_vcov(cal, ids)
  # the independence matrix reproduces the per-group-list results exactly
  le_ind <- linking_error(link, method = "jackknife_bc",
                          vcov = sv$vcov_independence)
  expect_equal(le_ind$le, le_list$le, tolerance = 1e-10)
  # the dependent matrix passes through and yields finite errors
  le_dep <- linking_error(link, method = "jackknife_bc", vcov = sv$vcov)
  expect_true(all(is.finite(le_dep$le$se_mu[-1])))
  expect_true(all(is.finite(le_dep$le$se_sigma[-1])))
})

test_that("stacked_vcov validates its inputs", {
  set.seed(435)
  d <- sim_trend_data(n_groups = 2, N = 200, I = 6, dif = "none")
  d <- lapply(d, function(X) as.matrix(as.data.frame(X)))
  cal <- suppressWarnings(calibrate(d, model = "2PL", calibration = "separate",
                   engine = "em"))
  expect_error(stacked_vcov(cal, list(1:200)), "one identifier vector")
  expect_error(stacked_vcov(cal, list(1:199, 1:200)), "one identifier per row")
  expect_error(stacked_vcov(cal, list(c(1, 1, 3:200), 1:200)), "duplicated")
  cal$models <- NULL
  expect_error(stacked_vcov(cal, list(1:200, 1:200)), "keep_models")
})
