# testlet (cluster) units in the jackknife and approximate jackknife

ip_testlet <- function() {
  mu_t <- c(0, 0.3, 0.6)
  sigma_t <- c(1, 1.1, 1.2)
  b_ref <- seq(-1.5, 1.5, length.out = 10)
  a_ref <- rep(c(0.8, 1.2), 5)
  ip <- do.call(rbind, lapply(1:3, function(t) {
    data.frame(group = t, item = sprintf("I%02d", 1:10),
               a = a_ref * sigma_t[t],
               b = (b_ref - mu_t[t]) / sigma_t[t])
  }))
  # shared drift on the first testlet at group 2
  ip$b[ip$group == 2 & ip$item %in% c("I01", "I02")] <-
    ip$b[ip$group == 2 & ip$item %in% c("I01", "I02")] + 0.4
  ip
}

testlet_map <- c(I01 = "T1", I02 = "T1", I03 = "T2", I04 = "T2",
                 I05 = "T3", I06 = "T3", I07 = "T4", I08 = "T4")

test_that("singleton clusters reproduce the item jackknife", {
  lk <- link_groups(ip_testlet(), approach = "chain", method = "mgm")
  jk <- linking_error(lk, method = "jackknife")
  singles <- stats::setNames(paste0("S", 1:10), sprintf("I%02d", 1:10))
  jk_s <- linking_error(lk, method = "jackknife", cluster = singles)
  expect_equal(jk$le$le_mu, jk_s$le$le_mu, tolerance = 1e-12)
  expect_equal(jk$le$le_sigma, jk_s$le$le_sigma, tolerance = 1e-12)
})

test_that("testlet units change the jackknife and count as n_reps", {
  lk <- link_groups(ip_testlet(), approach = "chain", method = "mgm")
  jk_i <- linking_error(lk, method = "jackknife")
  jk_c <- linking_error(lk, method = "jackknife", cluster = testlet_map)
  # 4 testlets plus 2 unmapped singleton items
  expect_true(all(jk_c$le$n_reps == 6))
  expect_true(all(jk_i$le$n_reps == 10))
  expect_gt(abs(jk_c$le$le_mu[2] - jk_i$le$le_mu[2]), 1e-8)
  # the shared drift sits in one testlet, deleting it jointly sees
  # the full shift
  expect_gt(jk_c$le$le_mu[2], jk_i$le$le_mu[2])
})

test_that("cluster ajk tracks the cluster jackknife for mgm and mm", {
  for (m in c("mgm", "mm")) {
    lk <- link_groups(ip_testlet(), approach = "chain", method = m)
    jk <- linking_error(lk, method = "jackknife", cluster = testlet_map)
    aj <- linking_error(lk, method = "ajk", cluster = testlet_map)
    expect_true(all(aj$le$n_reps == 6))
    expect_equal(aj$le$le_mu[-1], jk$le$le_mu[-1], tolerance = 0.1)
  }
})

test_that("cluster ajk runs for joint PHL", {
  lk <- link_groups(ip_testlet(), approach = "joint", method = "haberman",
                    variant = "pairwise")
  aj <- linking_error(lk, method = "ajk", cluster = testlet_map)
  expect_true(all(aj$le$n_reps == 6))
  expect_true(all(is.finite(aj$le$le_mu[-1])))
})

test_that("cluster guards fire", {
  lk <- link_groups(ip_testlet(), approach = "chain", method = "mgm")
  expect_error(linking_error(lk, method = "sandwich_esw",
                             cluster = testlet_map), "Cluster")
  expect_error(linking_error(lk, method = "jackknife",
                             cluster = unname(testlet_map)), "named")
  lh <- link_groups(ip_testlet(), approach = "chain", method = "haebara")
  expect_error(linking_error(lh, method = "ajk", cluster = testlet_map),
               "jackknife")
})

test_that("cluster jackknife works for response function and GPCM links", {
  lh <- link_groups(ip_testlet(), approach = "chain", method = "haebara")
  jk <- linking_error(lh, method = "jackknife", cluster = testlet_map)
  expect_true(all(jk$le$n_reps == 6))
  ipg <- ip_testlet()
  ipg$tau1 <- -0.4
  ipg$tau2 <- 0.4
  lg <- link_groups(ipg, approach = "chain", method = "sl")
  jkg <- linking_error(lg, method = "jackknife", cluster = testlet_map)
  expect_true(all(is.finite(jkg$le$le_mu[-1])))
})
