test_that("em_theta_grid spans the range with the requested nodes", {
  g <- em_theta_grid()
  expect_length(g, 101L)
  expect_equal(range(g), c(-6, 6))
  g61 <- em_theta_grid(nodes = 61L, range = c(-5, 5))
  expect_length(g61, 61L)
  expect_equal(range(g61), c(-5, 5))
})

test_that("em_irf_matrix computes clamped 2PL probabilities", {
  theta <- c(-1, 0, 2)
  P <- em_irf_matrix(a = c(1, 2), nu = c(0, 1), theta)
  expect_equal(dim(P), c(2L, 3L))
  expect_equal(P[1, ], plogis(theta), tolerance = 1e-12)
  expect_equal(P[2, ], plogis(2 * theta - 1), tolerance = 1e-12)
  Pex <- em_irf_matrix(a = 10, nu = 10, theta = c(-6, 6))
  expect_true(all(Pex >= 1e-10 & Pex <= 1 - 1e-10))
})
