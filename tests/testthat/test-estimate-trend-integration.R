# tests/testthat/test-estimate-trend-integration.R

test_that("all valid no-DIF matrix cells recover the same trend", {
  skip_on_cran()
  set.seed(4501)
  d <- sim_trend_data(n_groups = 3, N = 2500, I = 20, dif = "none",
                      mu = c(0, 0.4, 0.8), sigma = c(1, 1.1, 1.2))
  truth <- c(0, 0.4, 0.8)
  # link is consulted only for separate; mm/mgm are chain-only, so the joint
  # separate cells use haberman. For non-separate cells link is ignored.
  cells <- list(
    list("separate",    "chain",            "mgm"),
    list("separate",    "joint",            "haberman"),
    list("separate",    "joint_restricted", "haberman"),
    list("concurrent",  "chain",            "haberman"),
    list("concurrent",  "joint",            "haberman"),
    list("concurrent",  "joint_restricted", "haberman"),
    list("fixed",       "chain",            "haberman"),
    list("regularized", "chain",            "haberman")
  )
  for (cell in cells) {
    tr <- estimate_trend(d, calibration = cell[[1]], approach = cell[[2]],
                         link = cell[[3]])
    expect_s3_class(tr, "irtlink_trend")
    expect_equal(tr$trend$mu, truth, tolerance = 0.15,
                 info = paste(cell[[1]], cell[[2]], sep = "/"))
  }
})
