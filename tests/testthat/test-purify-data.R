test_that("data-based iterative purification flags the drifted items", {
  skip_on_cran()
  set.seed(1201)
  d <- sim_trend_data(n_groups = 3, N = 2500, I = 20, dif = "unbalanced",
                      dif_pct = 0.2, dif_effect = 0.7, mu = c(0, 0.4, 0.8),
                      sigma = c(1, 1, 1))
  drift <- attr(d, "dif_items")
  cal <- calibrate(d)
  dif <- detect_dif(cal, rmsd_method = "data", cutoff = 0.05,
                    procedure = "iterative_forward")
  expect_s3_class(dif, "irtlink_dif")
  expect_gt(length(intersect(dif$flagged, drift)), 0)
  expect_length(intersect(dif$anchor, drift), 0)
  expect_true(length(dif$history) >= 1)
})

test_that("data-based purification beats the naive one-step under heavy drift", {
  skip_on_cran()
  set.seed(1202)
  d <- sim_trend_data(n_groups = 3, N = 2500, I = 20, dif = "unbalanced",
                      dif_pct = 0.3, dif_effect = 0.8, mu = c(0, 0.5, 1.0),
                      sigma = c(1, 1, 1))
  cal <- calibrate(d)
  dif <- detect_dif(cal, rmsd_method = "data",
                    procedure = "iterative_forward")
  tr_pure  <- link_chain(cal, method = "mgm", anchor = dif$anchor)$trend
  tr_naive <- link_chain(cal, method = "mgm")$trend
  err_pure  <- max(abs(tr_pure$mu  - c(0, 0.5, 1.0)))
  err_naive <- max(abs(tr_naive$mu - c(0, 0.5, 1.0)))
  expect_lt(err_pure, err_naive)
})
