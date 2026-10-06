# tests/testthat/test-dif-methods.R

make_dif_fixture <- function() {
  tab <- data.frame(item = sprintf("I%02d", 1:5),
                    w2 = c(0.02, 0.30, 0.01, 0.04, 0.02),
                    w3 = c(0.03, 0.28, 0.02, 0.05, 0.03),
                    max = c(0.03, 0.30, 0.02, 0.05, 0.03))
  structure(
    list(rmsd = tab, rmsd_parameter = tab, rmsd_data = NULL,
         flagged = "I02", anchor = c("I01", "I03", "I04", "I05"),
         flagged_by_cutoff = list("0.05" = "I02"),
         linkable = tab$item, history = NULL, converged = TRUE,
         method = "rmsd", rmsd_method = "parameter", cutoff_type = "fixed",
         cutoff = 0.05, tau = 2.7, min_anchor = 3, procedure = "one_step",
         call = NULL),
    class = "irtlink_dif"
  )
}

test_that("print shows method, flagged count and anchor size", {
  res <- make_dif_fixture()
  expect_output(print(res), "RMSD")
  expect_output(print(res), "Flagged: 1")
  expect_output(print(res), "Anchor: 4")
})

test_that("summary lists flagged items and their statistics", {
  res <- make_dif_fixture()
  expect_output(summary(res), "I02")
})

test_that("plot draws RMSD with a cutoff line without error", {
  res <- make_dif_fixture()
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  expect_no_error(plot(res))
})

test_that("coef returns the driver RMSD table", {
  res <- make_dif_fixture()
  expect_identical(coef(res), res$rmsd)
})

test_that("plot handles an empty RMSD table without error", {
  res <- make_dif_fixture()
  res$rmsd <- data.frame(item = character(0), max = numeric(0))
  res$flagged <- character(0)
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  expect_no_error(plot(res))
})

test_that("print renders cleanly when nothing is flagged", {
  res <- make_dif_fixture()
  res$flagged <- character(0)
  res$anchor <- res$linkable
  expect_output(print(res), "Flagged: 0")
})

test_that("summary shows the per-cutoff block for multiple cutoffs", {
  res <- make_dif_fixture()
  res$flagged_by_cutoff <- list("0.05" = "I02", "0.10" = character(0))
  expect_output(summary(res), "Flagging by cutoff")
})

make_lrt_fixture <- function() {
  lt <- data.frame(item = sprintf("I%02d", 1:5),
                   chisq = c(1.2, 18.5, 0.8, 2.1, 1.0),
                   df = rep(2, 5),
                   p = c(0.55, 1e-4, 0.67, 0.35, 0.61),
                   p_adj = c(1, 5e-4, 1, 1, 1))
  structure(
    list(rmsd = NULL, rmsd_parameter = NULL, rmsd_data = NULL, lrt = lt,
         flagged = "I02", anchor = c("I01", "I03", "I04", "I05"),
         flagged_by_cutoff = list("0.05" = "I02"),
         linkable = lt$item, history = NULL, converged = TRUE,
         method = "lrt", rmsd_method = NA_character_,
         cutoff_type = NA_character_, cutoff = NA_real_, tau = NA_real_,
         alpha = 0.05, alpha_adjust = "none", dif_type = "both",
         min_anchor = 3, procedure = "one_step", call = NULL),
    class = "irtlink_dif"
  )
}

test_that("print shows LRT, dif_type and the smallest p-values", {
  res <- make_lrt_fixture()
  expect_output(print(res), "LRT")
  expect_output(print(res), "dif_type = both")
  expect_output(print(res), "Flagged: 1")
})

test_that("summary and coef handle an LRT object", {
  res <- make_lrt_fixture()
  expect_output(summary(res), "I02")
  expect_identical(coef(res), res$lrt)
})

test_that("plot draws an LRT object with a threshold line without error", {
  res <- make_lrt_fixture()
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  expect_no_error(plot(res))
})

test_that("print shows Bonferroni when adjusted", {
  res <- make_lrt_fixture()
  res$alpha_adjust <- "bonferroni"
  expect_output(print(res), "Bonferroni")
})
