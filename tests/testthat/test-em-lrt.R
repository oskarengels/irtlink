em_lrt_data <- function(seed = 95, drift = 0.6) {
  set.seed(seed)
  J <- 12L; N <- 2000L
  a <- runif(J, 0.7, 1.6); b <- rnorm(J)
  items <- paste0("i", seq_len(J))
  th1 <- rnorm(N, 0, 1); th2 <- rnorm(N, 0.2, 1)
  eta1 <- sweep(outer(th1, a), 2, a * b, "-")
  eta2 <- sweep(outer(th2, a), 2, a * b, "-")
  eta2[, 3] <- eta2[, 3] - drift          # uniform drift on item 3
  X1 <- matrix(rbinom(N * J, 1L, plogis(eta1)), N, J)
  X2 <- matrix(rbinom(N * J, 1L, plogis(eta2)), N, J)
  colnames(X1) <- colnames(X2) <- items
  list(X1, X2)
}

test_that("em LRT flags the drifted item and spares clean ones", {
  skip_on_cran()
  d <- em_lrt_data()
  cal <- calibrate(d, engine = "em")
  lt <- lrt_item_deviances(cal, dif_type = "both")
  expect_true(all(lt$chisq >= 0))
  expect_lt(lt$p[lt$item == "i3"], 0.001)
  expect_gt(min(lt$p[lt$item != "i3"]), 0.01)
  expect_equal(lt$df, rep(2L, nrow(lt)))
})

test_that("em LRT dif_type variants have the right df", {
  skip_on_cran()
  d <- em_lrt_data(96)
  cal <- calibrate(d, engine = "em")
  ltu <- lrt_item_deviances(cal, dif_type = "uniform")
  ltn <- lrt_item_deviances(cal, dif_type = "nonuniform")
  expect_equal(ltu$df, rep(1L, nrow(ltu)))
  expect_equal(ltn$df, rep(1L, nrow(ltn)))
  expect_lt(ltu$p[ltu$item == "i3"], 0.001)   # uniform drift -> uniform test
})

test_that("detect_dif lrt works on an em calibration", {
  skip_on_cran()
  d <- em_lrt_data(97)
  cal <- calibrate(d, engine = "em")
  dif <- detect_dif(cal, method = "lrt", alpha = 0.01)
  expect_true("i3" %in% dif$flagged)
})
