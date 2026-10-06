test_that("all design: shared link items plus group-unique items", {
  set.seed(1)
  d <- sim_trend_data(n_groups = 3, N = 50, I = 20, overlap = 0.4,
                      design = "all")
  expect_length(d, 3)
  expect_true(all(vapply(d, nrow, integer(1)) == 50))
  expect_true(all(vapply(d, ncol, integer(1)) == 20))
  # 8 anchor items (round(20 * 0.4)) present in every group
  common_all <- Reduce(intersect, lapply(d, colnames))
  expect_length(common_all, 8)
  # unique items appear in exactly one group
  all_items <- unlist(lapply(d, colnames))
  uniq <- setdiff(unique(all_items), common_all)
  expect_true(all(table(all_items[all_items %in% uniq]) == 1))
})

test_that("responses are dichotomous and attributes carry truth", {
  set.seed(1)
  mu <- c(0, 0.3, 0.6); sigma <- c(1, 1.1, 1.2)
  d <- sim_trend_data(n_groups = 3, N = 30, I = 10, mu = mu, sigma = sigma)
  vals <- unlist(lapply(d, function(x) unlist(x)))
  expect_true(all(vals %in% c(0, 1)))
  expect_equal(attr(d, "true_mu"), mu)
  expect_equal(attr(d, "true_sigma"), sigma)
  expect_s3_class(attr(d, "true_ipars"), "data.frame")
})

test_that("defaults: mu in 0.3 steps, sigma constant 1", {
  set.seed(1)
  d <- sim_trend_data(n_groups = 4, N = 10, I = 10)
  expect_equal(attr(d, "true_mu"), c(0, 0.3, 0.6, 0.9))
  expect_equal(attr(d, "true_sigma"), rep(1, 4))
})

test_that("invalid inputs error clearly", {
  expect_error(sim_trend_data(n_groups = 1, N = 10), "at least two groups")
  expect_error(sim_trend_data(n_groups = 3, N = 10, I = 10, overlap = 0.05),
               "at least 2")
})

test_that("mu/sigma/overlap validation gives clear errors", {
  expect_error(sim_trend_data(n_groups = 3, N = 10, mu = c(0, 0.3)),
               "`mu` must have length")
  expect_error(sim_trend_data(n_groups = 3, N = 10, sigma = c(1, 1)),
               "`sigma` must have length")
  expect_error(sim_trend_data(n_groups = 2, N = 10, sigma = c(1, 0)),
               "must be positive")
  expect_error(sim_trend_data(n_groups = 2, N = 10, overlap = 1.1),
               "must be in \\(0, 1\\]")
})

test_that("successive design: neighbors share link blocks, non-neighbors none", {
  set.seed(1)
  d <- sim_trend_data(n_groups = 4, N = 20, I = 20, overlap = 0.3,
                      design = "successive")
  n_link <- round(20 * 0.3)  # 6
  for (t in 1:3) {
    expect_length(intersect(colnames(d[[t]]), colnames(d[[t + 1]])), n_link)
  }
  expect_length(intersect(colnames(d[[1]]), colnames(d[[3]])), 0)
  expect_length(intersect(colnames(d[[1]]), colnames(d[[4]])), 0)
  expect_true(all(vapply(d, ncol, integer(1)) == 20))
})

test_that("successive design rejects overlap that exceeds group capacity", {
  expect_error(
    sim_trend_data(n_groups = 3, N = 10, I = 10, overlap = 0.6,
                   design = "successive"),
    "2 \\* round\\(I \\* overlap\\)"
  )
})

test_that("successive design: n_groups=2 allows overlap > 0.5", {
  set.seed(1)
  d <- sim_trend_data(n_groups = 2, N = 10, I = 10, overlap = 0.7,
                      design = "successive")
  expect_length(intersect(colnames(d[[1]]), colnames(d[[2]])), 7L)
  expect_true(all(vapply(d, ncol, integer(1)) == 10L))
})

test_that("unbalanced DIF shifts a known subset of link item difficulties", {
  set.seed(1)
  d <- sim_trend_data(n_groups = 3, N = 50, I = 20, overlap = 0.4,
                      design = "all", dif = "unbalanced",
                      dif_pct = 0.5, dif_effect = 0.8)
  flagged <- attr(d, "dif_items")
  # round(8 anchors * 0.5) = 4 affected items, all named
  expect_length(flagged, 4)
  expect_true(all(flagged %in% Reduce(intersect, lapply(d, colnames))))
  # base params unchanged; drift lives in attr
  expect_s3_class(attr(d, "true_ipars"), "data.frame")
  expect_identical(attr(d, "dif"), "unbalanced")
})

test_that("balanced DIF affects items in both directions, random uses dif_sd", {
  set.seed(1)
  db <- sim_trend_data(n_groups = 3, N = 30, I = 20, dif = "balanced",
                       dif_pct = 0.5, dif_effect = 0.8)
  expect_length(attr(db, "dif_items"), 4)
  set.seed(1)
  dr <- sim_trend_data(n_groups = 3, N = 30, I = 20, dif = "random",
                       dif_pct = 0.5, dif_sd = 0.3)
  expect_length(attr(dr, "dif_items"), 4)
  expect_identical(attr(dr, "dif"), "random")
})

test_that("dif = none leaves dif_items empty and still works", {
  set.seed(1)
  d <- sim_trend_data(n_groups = 3, N = 30, I = 10)
  expect_length(attr(d, "dif_items"), 0)
  expect_identical(attr(d, "dif"), "none")
})

test_that("DIF requires at least one affected item", {
  expect_error(
    sim_trend_data(n_groups = 3, N = 30, I = 20, dif = "unbalanced",
                   dif_pct = 0.01),
    "at least one"
  )
})

test_that("dif_effects attribute records the true per-group shifts", {
  set.seed(1)
  d <- sim_trend_data(n_groups = 3, N = 30, I = 20, dif = "unbalanced",
                      dif_pct = 0.5, dif_effect = 0.8)
  dr <- attr(d, "dif_effects")
  expect_true(is.matrix(dr))
  expect_equal(ncol(dr), 3)
  # group 1 is the reference: no drift
  expect_true(all(dr[, 1] == 0))
  # affected items carry +0.8 at groups 2 and 3
  drifted <- attr(d, "dif_items")
  expect_true(all(dr[drifted, 2] == 0.8))
  expect_true(all(dr[drifted, 3] == 0.8))
  # untouched items stay 0
  expect_true(all(dr[setdiff(rownames(dr), drifted), ] == 0))
})

test_that("balanced DIF warns when the affected count is odd", {
  expect_warning(
    sim_trend_data(n_groups = 2, N = 20, I = 10, overlap = 0.6,
                   dif = "balanced", dif_pct = 0.5, dif_effect = 0.8),
    "odd number"
  )
})

test_that("dif_effects attribute records the per-group drift magnitudes", {
  set.seed(1)
  d <- sim_trend_data(n_groups = 3, N = 30, I = 20, dif = "random",
                      dif_pct = 0.5, dif_sd = 0.3)
  dr <- attr(d, "dif_effects")
  expect_true(is.matrix(dr))
  expect_equal(ncol(dr), 3)
  # group 1 is the reference: no drift
  expect_true(all(dr[, 1] == 0))
  # drifted items have nonzero drift at >= one later group; others are all 0
  drifted <- attr(d, "dif_items")
  expect_true(all(rowSums(abs(dr[drifted, , drop = FALSE])) > 0))
  others <- setdiff(rownames(dr), drifted)
  expect_true(all(dr[others, ] == 0))
})

test_that("balanced DIF warns on an odd number of affected items", {
  expect_warning(
    sim_trend_data(n_groups = 2, N = 20, I = 10, dif = "balanced",
                   dif_pct = 0.75),
    "odd number"
  )
})
