# tests/testthat/test-detect-dif-integration.R

test_that("LRT and data-RMSD flag overlapping items under clear drift", {
  skip_on_cran()
  set.seed(1601)
  d <- sim_trend_data(n_groups = 2, N = 3000, I = 18, dif = "unbalanced",
                      dif_pct = 0.22, dif_effect = 0.9, mu = c(0, 0.4),
                      sigma = c(1, 1))
  drift <- attr(d, "dif_items")
  cal <- calibrate(d)
  lrt  <- detect_dif(cal, method = "lrt")
  data <- detect_dif(cal, rmsd_method = "data", cutoff = 0.05)
  # both recover a majority of the injected drift
  expect_gt(length(intersect(lrt$flagged,  drift)) / length(drift), 0.5)
  expect_gt(length(intersect(data$flagged, drift)) / length(drift), 0.5)
})

test_that("data-based iterative purification yields a near-true trend", {
  skip_on_cran()
  set.seed(1602)
  d <- sim_trend_data(n_groups = 3, N = 2500, I = 20, dif = "unbalanced",
                      dif_pct = 0.25, dif_effect = 0.8, mu = c(0, 0.4, 0.8),
                      sigma = c(1, 1, 1))
  cal <- calibrate(d)
  dif <- detect_dif(cal, rmsd_method = "data",
                    procedure = "iterative_forward")
  tr <- link_chain(cal, method = "mgm", anchor = dif$anchor)$trend
  expect_equal(tr$mu, c(0, 0.4, 0.8), tolerance = 0.15)
})
