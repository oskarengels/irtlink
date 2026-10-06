test_that("estimate_trend separate+purify beats the contaminated link", {
  skip_on_cran()
  set.seed(4301)
  d <- sim_trend_data(n_groups = 3, N = 2500, I = 20, dif = "unbalanced",
                      dif_pct = 0.25, dif_effect = 0.8, mu = c(0, 0.5, 1.0),
                      sigma = c(1, 1, 1))
  pure  <- estimate_trend(d, "separate", "chain", link = "mgm",
                          dif = "purify", procedure = "iterative_forward")
  naive <- estimate_trend(d, "separate", "chain", link = "mgm", dif = "none")
  expect_false(is.null(pure$dif))
  expect_lt(max(abs(pure$trend$mu  - c(0, 0.5, 1.0))),
            max(abs(naive$trend$mu - c(0, 0.5, 1.0))))
})

test_that("estimate_trend concurrent+joint+purify recovers under drift", {
  skip_on_cran()
  set.seed(4302)
  d <- sim_trend_data(n_groups = 3, N = 2500, I = 20, dif = "unbalanced",
                      dif_pct = 0.2, dif_effect = 0.7, mu = c(0, 0.4, 0.8),
                      sigma = c(1, 1, 1))
  tr <- estimate_trend(d, "concurrent", "joint", dif = "purify")
  expect_false(is.null(tr$dif))
  expect_equal(tr$trend$mu, c(0, 0.4, 0.8), tolerance = 0.15)
})

test_that("concurrent chain/restricted + purify is rejected", {
  d <- sim_trend_data(n_groups = 2, N = 300, I = 10, dif = "none")
  expect_error(estimate_trend(d, "concurrent", "chain", dif = "purify"),
               "joint")
})
