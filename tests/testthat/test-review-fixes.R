# regression tests for the external code review findings

test_that("person_overlap = 1 keeps the panel identifiers intact", {
  set.seed(3)
  d <- sim_trend_data(n_groups = 3, N = 50, I = 8, overlap = 1,
                      person_overlap = 1)
  expect_identical(rownames(d[[1]]), rownames(d[[2]]))
  expect_identical(rownames(d[[2]]), rownames(d[[3]]))
  expect_identical(vapply(d, nrow, integer(1)), rep(50L, 3))
})

test_that("data-driven flagging works with a single table row", {
  tab <- data.frame(item = "I01", g1 = 0.02, g2 = 0.08, max = 0.08)
  out <- irtlink:::dif_flag_mad(tab, tau = 2.7, min_anchor = 0)
  expect_type(out, "character")
})

test_that("bootstrap draws reuse returns the overlap matrix", {
  set.seed(4)
  d <- sim_trend_data(n_groups = 2, N = 200, I = 6, overlap = 1,
                      person_overlap = 0.5)
  ids <- attr(d, "person_ids")
  bv <- dependent_vcov(d, ids, method = "bootstrap", B = 4, seed = 1)
  bv2 <- dependent_vcov(d, ids, method = "bootstrap", draws = bv$draws)
  expect_false(is.null(bv2$overlap))
  expect_identical(bv2$overlap, bv$overlap)
})

test_that("LRT accepts control overrides through the dots", {
  set.seed(6)
  d <- sim_trend_data(n_groups = 2, N = 150, I = 6, overlap = 1)
  cal <- calibrate(d)
  res <- detect_dif(cal, method = "lrt",
                    control = list(conv = 1e-5, maxit = 300L))
  expect_s3_class(res, "irtlink_dif")
})

test_that("a failing delete-one refit is discarded, not fatal", {
  ip <- data.frame(
    group = c(1, 1, 1, 2, 2, 2, 3, 3),
    item = c("A", "B", "C", "A", "B", "C", "A", "D"),
    a = c(1, 1.2, 0.9, 1.1, 1.3, 1, 1.2, 0.8),
    b = c(-0.5, 0.2, 0.6, -0.7, 0, 0.4, -0.9, 0.3)
  )
  lk <- link_groups(ip, approach = "chain", method = "mgm")
  expect_warning(le <- linking_error(lk, method = "jackknife"),
                 "discarded")
  expect_true(is.finite(le$le$le_mu[2]))
})

test_that("plotting a joint link with indirectly connected groups works", {
  ip <- data.frame(
    group = c(1, 1, 1, 1, 2, 2, 3, 3),
    item = c("A", "B", "C", "D", "A", "B", "C", "D"),
    a = c(1, 1.1, 0.9, 1.2, 1, 1.15, 0.95, 1.25),
    b = c(-0.5, 0.2, 0.6, -0.2, -0.6, 0.1, 0.3, -0.4)
  )
  lk <- link_groups(ip, approach = "joint", method = "haberman")
  grDevices::pdf(NULL)
  expect_no_error(plot(lk))
  grDevices::dev.off()
})

test_that("the LRT plot handles item-specific degrees of freedom", {
  d_obj <- structure(
    list(method = "lrt",
         lrt = data.frame(item = c("A", "B", "C"), chisq = c(5, 3, 1),
                          df = c(4, 2, 2), p = c(0.1, 0.2, 0.5)),
         flagged = character(0), alpha = 0.05),
    class = "irtlink_dif")
  grDevices::pdf(NULL)
  expect_no_error(plot(d_obj))
  grDevices::dev.off()
})
