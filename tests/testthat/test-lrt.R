test_that("lrt_item_deviances separates drifted from anchor items (dif_type both)", {
  skip_on_cran()
  set.seed(1301)
  d <- sim_trend_data(n_groups = 2, N = 2500, I = 16, dif = "unbalanced",
                      dif_pct = 0.25, dif_effect = 0.8, mu = c(0, 0.4),
                      sigma = c(1, 1))
  drift  <- attr(d, "dif_items")
  cal <- calibrate(d)
  lt <- lrt_item_deviances(cal, dif_type = "both")
  expect_true(all(c("item", "chisq", "df", "p") %in% names(lt)))
  expect_true(all(lt$df == 2))            # both, 2 groups -> df = 2
  p_drift  <- lt$p[lt$item %in% drift]
  p_anchor <- lt$p[!lt$item %in% drift]
  expect_lt(stats::median(p_drift), stats::median(p_anchor))
})

test_that("dif_type controls the degrees of freedom", {
  skip_on_cran()
  set.seed(1302)
  d <- sim_trend_data(n_groups = 2, N = 1500, I = 14, dif = "unbalanced",
                      dif_pct = 0.25, dif_effect = 0.7, mu = c(0, 0.4),
                      sigma = c(1, 1))
  cal <- calibrate(d)
  expect_true(all(lrt_item_deviances(cal, dif_type = "uniform")$df == 1))
  expect_true(all(lrt_item_deviances(cal, dif_type = "nonuniform")$df == 1))
})

test_that("detect_dif(method = 'lrt') flags drifted items; bonferroni is stricter", {
  skip_on_cran()
  set.seed(1401)
  d <- sim_trend_data(n_groups = 2, N = 2500, I = 16, dif = "unbalanced",
                      dif_pct = 0.25, dif_effect = 0.8, mu = c(0, 0.4),
                      sigma = c(1, 1))
  drift <- attr(d, "dif_items")
  cal <- calibrate(d)
  dif <- detect_dif(cal, method = "lrt", alpha = 0.05)
  expect_s3_class(dif, "irtlink_dif")
  expect_identical(dif$method, "lrt")
  expect_gt(length(intersect(dif$flagged, drift)), 0)
  expect_length(intersect(dif$anchor, drift), 0)
  expect_false(is.null(dif$lrt))                     # stats table stored
  expect_identical(dif$dif_type, "both")             # default
  dif_b <- detect_dif(cal, method = "lrt", alpha = 0.05,
                      alpha_adjust = "bonferroni")
  expect_true(length(dif_b$flagged) <= length(dif$flagged))   # stricter
})

test_that("detect_dif LRT honours dif_type = uniform", {
  skip_on_cran()
  set.seed(1402)
  d <- sim_trend_data(n_groups = 2, N = 2500, I = 16, dif = "unbalanced",
                      dif_pct = 0.25, dif_effect = 0.8, mu = c(0, 0.4),
                      sigma = c(1, 1))
  cal <- calibrate(d)
  dif <- detect_dif(cal, method = "lrt", dif_type = "uniform")
  expect_identical(dif$dif_type, "uniform")
  expect_true(all(dif$lrt$df == 1))
})

test_that("LRT handles a successive design with per-item degrees of freedom", {
  skip_on_cran()
  set.seed(1303)
  d <- sim_trend_data(n_groups = 3, N = 2000, I = 18, design = "successive",
                      dif = "none", mu = c(0, 0.3, 0.6), sigma = c(1, 1, 1))
  cal <- calibrate(d)
  lt <- expect_no_error(lrt_item_deviances(cal, dif_type = "both"))
  expect_gt(nrow(lt), 0)
  # df must reflect each item's own group count: 2*(n_present - 1) for "both"
  presence <- table(unlist(lapply(d, function(w) colnames(as.data.frame(w)))))
  np <- as.integer(presence[lt$item])
  expect_equal(lt$df, 2L * (np - 1L))
})
