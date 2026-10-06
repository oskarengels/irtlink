# tests/testthat/test-rmsd-data.R
# Data-based RMSD on the internal EM multiple group fit; keep N/I modest.

test_that("rmsd_data returns per-group RMSD plus bias-corrected column", {
  set.seed(42)
  d <- sim_trend_data(n_groups = 2, N = 800, I = 12, overlap = 0.5,
                      dif = "none")
  cal <- calibrate(d, model = "2PL")
  tab <- rmsd_data(cal)
  expect_true(all(c("item", "max", "max_bc") %in% names(tab)))
  # no drift -> all RMSD small
  expect_lt(max(tab$max), 0.06)
})

test_that("rmsd_data flags a drifted anchor item", {
  set.seed(42)
  d <- sim_trend_data(n_groups = 2, N = 1500, I = 12, overlap = 0.5,
                      dif = "unbalanced", dif_pct = 0.34, dif_effect = 1.0)
  drifted <- attr(d, "dif_items")
  cal <- calibrate(d, model = "2PL")
  tab <- rmsd_data(cal)
  worst <- tab$item[which.max(tab$max)]
  expect_true(worst %in% drifted)
})

test_that("rmsd_data errors when models were dropped", {
  set.seed(42)
  d <- sim_trend_data(n_groups = 2, N = 300, I = 10, overlap = 0.5)
  cal <- calibrate(d, keep_models = FALSE)
  expect_error(rmsd_data(cal), "keep_models")
})

test_that("rmsd_data handles 3 groups with w3 columns", {
  set.seed(7)
  d <- sim_trend_data(n_groups = 3, N = 600, I = 12, overlap = 0.5,
                      dif = "none")
  cal <- calibrate(d, model = "2PL")
  tab <- rmsd_data(cal)
  expect_true(all(c("g1", "g2", "g3", "g3_bc", "max", "max_bc") %in%
                    names(tab)))
})
