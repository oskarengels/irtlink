# tests/testthat/test-dif-flag.R

test_that("dif_flag_fixed flags items above the cutoff", {
  tab <- data.frame(item = c("A", "B", "C", "D", "E"),
                    max = c(0.02, 0.10, 0.01, 0.07, 0.03))
  expect_setequal(dif_flag_fixed(tab, cutoff = 0.05, min_anchor = 2),
                  c("B", "D"))
})

test_that("tie-breaking keeps anchors >= min_anchor by dropping lowest flags", {
  tab <- data.frame(item = c("A", "B", "C", "D", "E"),
                    max = c(0.20, 0.18, 0.16, 0.14, 0.01))
  # 4 items above 0.05 but min_anchor = 3 allows only 5 - 3 = 2 flags
  flagged <- dif_flag_fixed(tab, cutoff = 0.05, min_anchor = 3)
  expect_setequal(flagged, c("A", "B"))
})

test_that("no flags when nothing exceeds the cutoff", {
  tab <- data.frame(item = c("A", "B"), max = c(0.01, 0.02))
  expect_length(dif_flag_fixed(tab, cutoff = 0.05, min_anchor = 1), 0)
})

test_that("dif_tiebreak returns at most total - min_anchor items", {
  stat <- c(A = 0.2, B = 0.18, C = 0.16, D = 0.14)
  keep <- dif_tiebreak(stat, above = names(stat), min_anchor = 3, total = 4)
  expect_length(keep, 1)
  expect_identical(keep, "A")
})

test_that("dif_flag_mad flags robust-z outliers per group", {
  # group columns g1, g2; one clear outlier in g2
  tab <- data.frame(
    item = sprintf("I%02d", 1:8),
    g1 = c(0.02, 0.03, 0.02, 0.025, 0.03, 0.02, 0.028, 0.022),
    g2 = c(0.02, 0.03, 0.02, 0.025, 0.30, 0.02, 0.028, 0.022)
  )
  flagged <- dif_flag_mad(tab, tau = 2.7, min_anchor = 3,
                          group_cols = c("g1", "g2"))
  expect_identical(flagged, "I05")
})

test_that("dif_flag_mad returns nothing when the spread is uniform", {
  tab <- data.frame(item = c("A", "B", "C", "D"),
                    g1 = c(0.02, 0.021, 0.019, 0.020),
                    g2 = c(0.02, 0.021, 0.019, 0.020))
  expect_length(dif_flag_mad(tab, tau = 2.7, min_anchor = 1,
                             group_cols = c("g1", "g2")), 0)
})

test_that("dif_flag_mad zero-MAD guard yields no flags for identical values", {
  # all-identical group values -> MAD = 0 -> robust_z returns all zeros
  tab <- data.frame(item = c("A", "B", "C", "D"),
                    g1 = rep(0.05, 4), g2 = rep(0.05, 4))
  expect_length(dif_flag_mad(tab, tau = 2.7, min_anchor = 1,
                             group_cols = c("g1", "g2")), 0)
})
