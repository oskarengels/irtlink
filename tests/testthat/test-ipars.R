test_that("as_ipars validates required columns", {
  expect_error(as_ipars(data.frame(group = 1, item = "I1")),
               "Missing columns.*a, b")
  expect_error(as_ipars(list(1, 2)), "irtlink_calib.*data.frame")
})

test_that("as_ipars requires consecutive groups starting at 1", {
  bad <- data.frame(group = c(1, 3), item = c("I1", "I1"),
                    a = c(1, 1), b = c(0, 0))
  expect_error(as_ipars(bad), "consecutive integers")
})

test_that("as_ipars returns normalized data.frame", {
  x <- data.frame(group = c(2, 1), item = c("I1", "I1"),
                  a = c(1.2, 1.0), b = c(0.5, 0.0),
                  extra = c("x", "y"))
  res <- as_ipars(x)
  expect_identical(names(res), c("group", "item", "a", "b", "c"))
  expect_identical(res$group, c(1L, 2L))
})

test_that("ipars_pair returns common items in wide format", {
  ip <- data.frame(
    group = c(1, 1, 1, 2, 2, 2),
    item = c("I1", "I2", "I3", "I2", "I3", "I4"),
    a = 1:6 / 2, b = seq(-1, 1.5, by = 0.5)
  )
  pair <- ipars_pair(as_ipars(ip), 1, 2)
  expect_identical(pair$item, c("I2", "I3"))
  expect_identical(names(pair), c("item", "a_from", "b_from", "a_to", "b_to", "c_from", "c_to"))
  expect_equal(pair$a_from, c(1.0, 1.5))
  expect_equal(pair$a_to, c(2.0, 2.5))
})

test_that("ipars_pair errors when no common items exist", {
  ip <- data.frame(group = c(1, 2), item = c("I1", "I2"),
                   a = c(1, 1), b = c(0, 0))
  expect_error(ipars_pair(as_ipars(ip), 1, 2), "No common items")
})

test_that("as_ipars rejects NA in group, a, or b", {
  base <- data.frame(group = c(1, 2), item = c("I1", "I1"),
                     a = c(1, 1), b = c(0, 0))
  na_a <- base; na_a$a[1] <- NA
  expect_error(as_ipars(na_a), "must not contain NA")
  na_b <- base; na_b$b[2] <- NA
  expect_error(as_ipars(na_b), "must not contain NA")
  na_group <- base; na_group$group[1] <- NA
  expect_error(as_ipars(na_group), "must not contain NA")
})

test_that("as_ipars rejects duplicate (group, item) combinations", {
  dup <- data.frame(group = c(1, 1, 2), item = c("I1", "I1", "I1"),
                    a = c(1, 1.1, 1), b = c(0, 0.1, 0))
  expect_error(as_ipars(dup), "Duplicate \\(group, item\\)")
})

test_that("as_ipars converts factor items to character labels", {
  x <- data.frame(group = c(1, 2), item = factor(c("I1", "I1")),
                  a = c(1, 1), b = c(0, 0))
  expect_identical(as_ipars(x)$item, c("I1", "I1"))
})
