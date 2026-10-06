# tests/testthat/test-bootstrap-vcov.R

boot_test_data <- function() {
  set.seed(42)
  d <- sim_trend_data(n_groups = 2, N = 120, I = 8, overlap = 1,
                      dif = "none", mu = c(0, 0.3), sigma = c(1, 1.1))
  ids <- list(as.character(1:120), as.character(61:180))
  list(d = d, ids = ids)
}

# Reference implementation of the resampling scheme of the KESS person
# bootstrap script: set.seed(seed + b), draw the union sample with
# replacement, keep the rows of the drawn persons per group, drop items
# without variance, calibrate separately.
boot_reference_draws <- function(resp_list, ids, B, seed) {
  resp_list <- lapply(resp_list, as.matrix)
  ids <- lapply(ids, as.character)
  union_ids <- unique(unlist(ids))
  N_union <- length(union_ids)
  row_of <- lapply(ids, function(v) stats::setNames(seq_along(v), v))
  out <- lapply(seq_len(B), function(b) {
    set.seed(seed + b)
    draw <- sample(union_ids, N_union, replace = TRUE)
    resp_b <- lapply(seq_along(resp_list), function(k) {
      rows <- row_of[[k]][draw]
      rows <- rows[!is.na(rows)]
      X <- resp_list[[k]][rows, , drop = FALSE]
      rownames(X) <- NULL
      v <- apply(X, 2, function(col) stats::var(col, na.rm = TRUE))
      drop <- is.na(v) | v == 0
      list(X = X[, !drop, drop = FALSE], dropped = colnames(X)[drop])
    })
    dropped <- unlist(lapply(resp_b, `[[`, "dropped"))
    cal_b <- suppressWarnings(calibrate(lapply(resp_b, `[[`, "X"),
                                        model = "2PL",
                                        calibration = "separate",
                                        engine = "em"))
    ip <- cal_b$ipars
    data.frame(rep = b, group = ip$group, item = ip$item,
               a = ip$a, b = ip$b,
               converged = paste(cal_b$converged, collapse = ","),
               dropped = paste(dropped, collapse = ";"),
               stringsAsFactors = FALSE)
  })
  do.call(rbind, out)
}

test_that("bootstrap_vcov reproduces the script resampling exactly", {
  td <- boot_test_data()
  bv <- bootstrap_vcov(td$d, td$ids, B = 3, seed = 20260919)
  ref <- boot_reference_draws(td$d, td$ids, B = 3, seed = 20260919)
  expect_equal(bv$draws$a, ref$a)
  expect_equal(bv$draws$b, ref$b)
  expect_identical(bv$draws$item, ref$item)
  expect_identical(bv$draws$converged, ref$converged)
  expect_identical(bv$draws$dropped, ref$dropped)
})

test_that("bootstrap_vcov returns a usable covariance", {
  td <- boot_test_data()
  bv <- bootstrap_vcov(td$d, td$ids, B = 4, seed = 1)
  expect_s3_class(bv, "irtlink_bootstrap_vcov")
  expect_true(isSymmetric(bv$vcov))
  first_item <- colnames(as.matrix(td$d[[1]]))[1]
  expect_identical(rownames(bv$vcov)[1:2],
                   paste0("g1:", first_item, ":", c("a", "b")))
  expect_identical(bv$overlap["g1", "g2"], 60L)
  expect_output(print(bv), "Person-bootstrap covariance")

  cal <- calibrate(td$d, model = "2PL", calibration = "separate")
  lk <- link_chain(cal, method = "mgm")
  le <- linking_error(lk, method = "jackknife_bc", vcov = bv$vcov)$le
  expect_true(all(is.finite(le$se_mu)))
})

test_that("the parallel path reproduces the sequential draws", {
  skip_on_cran()
  ns_path <- getNamespaceInfo("irtlink", "path")
  if (file.exists(file.path(ns_path, "R", "bootstrap_vcov.R"))) {
    skip("parallel workers load the installed irtlink, not the dev copy")
  }
  td <- boot_test_data()
  bv_seq <- bootstrap_vcov(td$d, td$ids, B = 4, seed = 9)
  bv_par <- bootstrap_vcov(td$d, td$ids, B = 4, seed = 9, cores = 2)
  expect_equal(bv_par$draws, bv_seq$draws)
  expect_equal(bv_par$vcov, bv_seq$vcov)
})

test_that("bootstrap_vcov assembles the covariance from stored draws", {
  td <- boot_test_data()
  bv <- bootstrap_vcov(td$d, td$ids, B = 3, seed = 7)
  bv2 <- bootstrap_vcov(td$d, td$ids, draws = bv$draws)
  expect_equal(bv2$vcov, bv$vcov)
  expect_identical(bv2$B_used, bv$B_used)
  expect_identical(bv2$overlap, bv$overlap)
})

test_that("items without variance in a replicate are recorded as dropped", {
  set.seed(31)
  d <- sim_trend_data(n_groups = 2, N = 120, I = 8, overlap = 1,
                      dif = "none")
  # one item is solved by three persons, resamples missing them
  # drop the item
  d[[1]][[1]] <- c(1L, 1L, 1L, rep(0L, 117))
  d[[2]][[1]] <- c(1L, 1L, 1L, rep(0L, 117))
  ids <- list(as.character(1:120), as.character(61:180))
  bv <- bootstrap_vcov(d, ids, B = 8, seed = 2)
  expect_true(any(bv$draws$dropped != ""))
  expect_true(length(bv$excluded) > 0)
})

test_that("bootstrap_vcov excludes nonconverged and incomplete replicates", {
  set.seed(5)
  dr <- expand.grid(rep = 1:4, group = 1:2, item = c("I01", "I02"),
                    stringsAsFactors = FALSE)
  dr$a <- exp(rnorm(nrow(dr), sd = 0.1))
  dr$b <- rnorm(nrow(dr))
  dr$converged <- "TRUE,TRUE"
  dr$dropped <- ""
  dr$converged[dr$rep == 2] <- "TRUE,FALSE"
  dr$dropped[dr$rep == 4] <- "I01"
  bv2 <- bootstrap_vcov(ids = list("p1", "p2"), draws = dr)
  expect_equal(bv2$excluded, c(2, 4))
  expect_identical(bv2$B_used, 2L)
})

test_that("bootstrap_vcov validates its inputs", {
  td <- boot_test_data()
  expect_error(bootstrap_vcov(td$d, td$ids[1]), "one identifier vector")
  expect_error(bootstrap_vcov(td$d, td$ids, B = 1), "at least 2")
  expect_error(bootstrap_vcov(td$d, td$ids, B = 3, cores = 2),
               "need `seed`")
  bad_ids <- td$ids
  bad_ids[[1]][2] <- bad_ids[[1]][1]
  expect_error(bootstrap_vcov(td$d, bad_ids, B = 3, seed = 1),
               "duplicated")
  expect_error(bootstrap_vcov(ids = list("p1", "p2"),
                              draws = data.frame(rep = 1)),
               "columns")
})
