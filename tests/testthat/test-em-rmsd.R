# The exact comparison against the published reference implementation
# lives in the development oracle tests (tests_dev/).

test_that("em_rmsd_data flags a drifted item on an em calibration", {
  skip_on_cran()
  set.seed(102)
  J <- 12L; N <- 2000L
  a <- runif(J, 0.7, 1.6); b <- rnorm(J)
  items <- paste0("i", seq_len(J))
  th1 <- rnorm(N, 0, 1); th2 <- rnorm(N, 0.2, 1)
  eta1 <- sweep(outer(th1, a), 2, a * b, "-")
  eta2 <- sweep(outer(th2, a), 2, a * b, "-")
  eta2[, 3] <- eta2[, 3] - 0.8            # uniform drift on item 3
  X1 <- matrix(rbinom(N * J, 1L, plogis(eta1)), N, J)
  X2 <- matrix(rbinom(N * J, 1L, plogis(eta2)), N, J)
  colnames(X1) <- colnames(X2) <- items
  cal <- calibrate(list(X1, X2), engine = "em")
  tab <- em_rmsd_data(cal)
  expect_named(tab, c("item", "g1", "g1_bc", "g2", "g2_bc", "max",
                      "max_bc"))
  expect_gt(tab$max[tab$item == "i3"], 0.04)
  expect_gt(tab$max[tab$item == "i3"], max(tab$max[tab$item != "i3"]))
})

test_that("rmsd_data dispatches on the calibration engine", {
  skip_on_cran()
  set.seed(103)
  d <- sim_trend_data(n_groups = 2, N = 500, I = 12, dif = "none")
  cal <- calibrate(d, engine = "em")
  tab <- rmsd_data(cal)
  expect_true(all(c("max", "max_bc") %in% names(tab)))
})

test_that("em_rmsd_data_split drops the split items", {
  skip_on_cran()
  set.seed(104)
  d <- sim_trend_data(n_groups = 2, N = 500, I = 12, dif = "none")
  dl <- lapply(d, function(x) as.matrix(as.data.frame(x)))
  it <- colnames(dl[[1]])[1]
  tab <- em_rmsd_data_split(dl, split_items = it)
  expect_false(it %in% tab$item)
  expect_true(all(c("max", "max_bc") %in% names(tab)))
})


test_that("the bias-corrected RMSD removes most finite sample bias", {
  skip_on_cran()
  set.seed(99)
  d <- sim_trend_data(n_groups = 2, N = 500, I = 12, overlap = 1,
                      dif = "none")
  fit <- irtlink:::em_fit_multi(lapply(d, as.matrix))
  tab <- irtlink:::em_rmsd_table(fit)
  expect_lt(median(tab$max_bc), 0.6 * median(tab$max))
  set.seed(7)
  d2 <- sim_trend_data(n_groups = 2, N = 1500, I = 12, overlap = 1,
                       dif = "unbalanced", dif_pct = 0.1,
                       dif_effect = 0.8)
  fit2 <- irtlink:::em_fit_multi(lapply(d2, as.matrix))
  tab2 <- irtlink:::em_rmsd_table(fit2)
  difit <- attr(d2, "dif_items")
  expect_gt(tab2$max_bc[tab2$item == difit], 0.04)
})