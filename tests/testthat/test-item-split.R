test_that("split_items_resp creates per-group columns for split items only", {
  X <- matrix(c(1,0,1, 0,1,1, 1,1,0, 0,0,1), nrow = 4, byrow = TRUE,
              dimnames = list(NULL, c("I1","I2","I3")))
  group <- c(1L, 1L, 2L, 2L)
  sp <- split_items_resp(X, group, split_items = "I2", n_groups = 2)
  expect_true(all(c("I1","I3","I2__G1","I2__G2") %in% colnames(sp)))
  expect_false("I2" %in% colnames(sp))
  expect_equal(sp[group == 2, "I2__G1"], c(NA_real_, NA_real_))
  expect_equal(sp[group == 1, "I2__G2"], c(NA_real_, NA_real_))
  expect_equal(sp[group == 1, "I2__G1"], X[group == 1, "I2"])
})

test_that("split_items_resp with no split items returns the matrix unchanged", {
  X <- matrix(c(1,0,1,0), nrow = 2, dimnames = list(NULL, c("A","B")))
  sp <- split_items_resp(X, group = c(1L,2L), split_items = character(0),
                         n_groups = 2)
  expect_setequal(colnames(sp), c("A","B"))
})
