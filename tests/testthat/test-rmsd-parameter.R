# tests/testthat/test-rmsd-parameter.R

test_that("rmsd_parameter is ~0 for drift-free identified parameters", {
  ip <- as_ipars(make_drift_ipars(c(0, 0.3, 0.6), c(1, 1.1, 1.2)))
  tab <- rmsd_parameter(ip, link_method = "mgm")
  expect_true(all(c("item", "max") %in% names(tab)))
  expect_lt(max(tab$max), 1e-6)
})

test_that("rmsd_parameter flags the drifted item with the largest value", {
  ip <- as_ipars(make_drift_ipars(c(0, 0.3, 0.6), c(1, 1.1, 1.2),
                                  drift_items = "I05", drift = 1.0))
  tab <- rmsd_parameter(ip, link_method = "mgm")
  worst <- tab$item[which.max(tab$max)]
  expect_identical(worst, "I05")
  expect_gt(max(tab$max), 0.05)
})

test_that("rmsd_parameter respects an explicit anchor for linking", {
  ip <- as_ipars(make_drift_ipars(c(0, 0.3), c(1, 1.1),
                                  drift_items = "I01", drift = 1.0))
  tab <- rmsd_parameter(ip, link_method = "mgm",
                        anchor = sprintf("I%02d", 3:20))
  expect_identical(tab$item[which.max(tab$max)], "I01")
})

test_that("rmsd_parameter drops items administered in a single group", {
  ip2 <- as_ipars(make_drift_ipars(c(0, 0.3, 0.6), c(1, 1.1, 1.2), I = 10))
  unique_item <- data.frame(group = 1L, item = "U01", a = 1, b = 0, c = 0)
  ip2 <- as_ipars(rbind(ip2[names(unique_item)], unique_item))
  tab <- rmsd_parameter(ip2, link_method = "mgm")
  expect_false("U01" %in% tab$item)
  expect_true(all(tab$item %in% sprintf("I%02d", 1:10)))
  expect_true(nrow(tab) == 10)
})
