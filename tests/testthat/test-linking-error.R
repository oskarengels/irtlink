test_that("jackknife linking error is zero for identified drift-free chain links", {
  ip <- make_identified_ipars(c(0, 0.3, 0.6), c(1, 1.1, 1.2), I = 12)
  lk <- link_chain(ip, method = "mgm")
  le <- linking_error(lk, method = "jackknife")

  expect_s3_class(le, "irtlink_link")
  expect_s3_class(le$le, "data.frame")
  expect_named(le$le, c("group", "le_mu", "le_sigma", "method", "n_reps"))
  expect_lt(max(abs(le$le$le_mu), na.rm = TRUE), 1e-10)
  expect_lt(max(abs(le$le$le_sigma), na.rm = TRUE), 1e-10)
})

test_that("jackknife linking error reacts to item-level drift", {
  ip <- make_identified_ipars(c(0, 0.3), c(1, 1.1), I = 12)
  ip$b[ip$group == 2 & ip$item == "I01"] <- ip$b[ip$group == 2 & ip$item == "I01"] + 0.5
  lk <- link_chain(ip, method = "mgm")
  le <- linking_error(lk)

  expect_gt(le$le$le_mu[le$le$group == 2], 0)
  expect_true(all(le$le$n_reps == 12))
})

test_that("jackknife linking error works for joint PHL links", {
  ip <- make_identified_ipars(c(0, 0.3, 0.6), c(1, 1.1, 1.2), I = 12)
  lk <- link_joint(ip, method = "phl")
  le <- linking_error(lk)

  expect_lt(max(abs(le$le$le_mu), na.rm = TRUE), 1e-10)
  expect_lt(max(abs(le$le$le_sigma), na.rm = TRUE), 1e-10)
})

test_that("AJK linking error is zero for identified drift-free links", {
  ip <- make_identified_ipars(c(0, 0.3, 0.6), c(1, 1.1, 1.2), I = 12)

  for (link_method in c("mgm", "mm")) {
    lk <- link_chain(ip, method = link_method)
    le <- linking_error(lk, method = "ajk")

    expect_lt(max(abs(le$le$le_mu), na.rm = TRUE), 1e-8)
    expect_lt(max(abs(le$le$le_sigma), na.rm = TRUE), 1e-8)
    expect_identical(le$le$method, rep("ajk", nrow(le$le)))
  }

  for (iw in c("uniform", "inverse_admin")) {
    lk <- link_joint(ip, method = "phl", item_weights = iw)
    le <- linking_error(lk, method = "ajk")

    expect_lt(max(abs(le$le$le_mu), na.rm = TRUE), 1e-8)
    expect_lt(max(abs(le$le$le_sigma), na.rm = TRUE), 1e-8)
    expect_identical(le$le$method, rep("ajk", nrow(le$le)))
  }
})

test_that("AJK linking error approximates delete-one-item jackknife", {
  ip <- make_identified_ipars(c(0, 0.3), c(1, 1.1), I = 12)
  ip$b[ip$group == 2 & ip$item == "I01"] <- ip$b[ip$group == 2 & ip$item == "I01"] + 0.5
  lk <- link_chain(ip, method = "mgm")
  jk <- linking_error(lk, method = "jackknife")
  ajk <- linking_error(lk, method = "ajk")

  expect_gt(ajk$le$le_mu[ajk$le$group == 2], 0)
  expect_equal(ajk$le$le_mu, jk$le$le_mu, tolerance = 1e-10)
  expect_equal(ajk$le$le_sigma, jk$le$le_sigma, tolerance = 1e-10)
  expect_true(all(ajk$le$n_reps == 12))
})

test_that("analytic AJK derivatives match the numerical ones internally", {
  # `ajk` picks the analytic derivatives automatically; this verifies the
  # numerical estimating-equation derivatives against them.
  ip <- make_identified_ipars(c(0, 0.3, 0.6), c(1, 1.1, 1.2), I = 12)
  ip$b[ip$group == 2 & ip$item == "I01"] <- ip$b[ip$group == 2 & ip$item == "I01"] + 0.5
  ip$a[ip$group == 2 & ip$item == "I02"] <- ip$a[ip$group == 2 & ip$item == "I02"] * 1.2
  ip$b[ip$group == 3 & ip$item == "I03"] <- ip$b[ip$group == 3 & ip$item == "I03"] - 0.4
  ip$a[ip$group == 3 & ip$item == "I04"] <- ip$a[ip$group == 3 & ip$item == "I04"] * 0.85

  for (link_method in c("mgm", "mm")) {
    lk <- link_chain(ip, method = link_method)
    le_num <- le_ajk_chain(lk, method = "ajk", analytic = FALSE)
    le_an <- le_ajk_chain(lk, method = "ajk", analytic = TRUE)

    expect_equal(le_an$le_mu, le_num$le_mu, tolerance = 1e-8)
    expect_equal(le_an$le_sigma, le_num$le_sigma, tolerance = 1e-8)
  }

  for (iw in c("uniform", "inverse_admin")) {
    lk <- link_joint(ip, method = "phl", item_weights = iw)
    le_num <- le_ajk_phl(lk, method = "ajk", analytic = FALSE)
    le_an <- le_ajk_phl(lk, method = "ajk", analytic = TRUE)

    expect_equal(le_an$le_mu, le_num$le_mu, tolerance = 1e-6)
    expect_equal(le_an$le_sigma, le_num$le_sigma, tolerance = 1e-6)
  }
})

test_that("jackknife and AJK use Robitzsch full-centered variants", {
  ip <- make_identified_ipars(c(0, 0.3), c(1, 1.1), I = 12)
  ip$b[ip$group == 2 & ip$item == "I01"] <- ip$b[ip$group == 2 & ip$item == "I01"] + 0.5
  lk <- link_chain(ip, method = "mgm")
  jk <- linking_error(lk, method = "jackknife")
  ajk <- linking_error(lk, method = "ajk")

  expect_gt(jk$le$le_mu[jk$le$group == 2], 0)
  expect_equal(ajk$le$le_mu, jk$le$le_mu, tolerance = 1e-10)
  expect_equal(ajk$le$le_sigma, jk$le$le_sigma, tolerance = 1e-10)
  expect_identical(jk$le$method, rep("jackknife", nrow(jk$le)))
  expect_identical(ajk$le$method, rep("ajk", nrow(ajk$le)))
})

test_that("bias-corrected jackknife variants report SE, LEbc, TE, and TEbc", {
  ip <- make_identified_ipars(c(0, 0.3), c(1, 1.1), I = 6)
  ip$b[ip$group == 2 & ip$item == "I01"] <- ip$b[ip$group == 2 & ip$item == "I01"] + 0.5
  lk <- link_chain(ip, method = "mgm")
  V <- matrix(0, nrow = 2 * nrow(lk$ipars), ncol = 2 * nrow(lk$ipars))

  for (method in c("jackknife_bc", "ajk_bc")) {
    le <- linking_error(lk, method = method, vcov = V)

    expect_named(le$le, c(
      "group", "le_mu", "le_sigma", "method", "n_reps",
      "se_mu", "se_sigma", "lebc_mu", "lebc_sigma",
      "te_mu", "te_sigma", "tebc_mu", "tebc_sigma"
    ))
    expect_identical(le$le$method, rep(method, nrow(le$le)))
    expect_equal(le$le$se_mu, rep(0, nrow(le$le)))
    expect_equal(le$le$se_sigma, rep(0, nrow(le$le)))
    expect_equal(le$le$lebc_mu, le$le$le_mu, tolerance = 1e-10)
    expect_equal(le$le$lebc_sigma, le$le$le_sigma, tolerance = 1e-10)
    expect_equal(le$le$te_mu, le$le$le_mu, tolerance = 1e-10)
    expect_equal(le$le$te_sigma, le$le$le_sigma, tolerance = 1e-10)
    expect_equal(le$le$tebc_mu, le$le$le_mu, tolerance = 1e-10)
    expect_equal(le$le$tebc_sigma, le$le$le_sigma, tolerance = 1e-10)
  }
})

test_that("bias-corrected AJK works for PHL links", {
  ip <- make_identified_ipars(c(0, 0.3, 0.6), c(1, 1.1, 1.2), I = 6)
  ip$b[ip$group == 2 & ip$item == "I01"] <- ip$b[ip$group == 2 & ip$item == "I01"] + 0.5
  ip$a[ip$group == 3 & ip$item == "I02"] <- ip$a[ip$group == 3 & ip$item == "I02"] * 1.2
  lk <- link_joint(ip, method = "phl", item_weights = "inverse_admin")
  V <- matrix(0, nrow = 2 * nrow(lk$ipars), ncol = 2 * nrow(lk$ipars))
  le <- linking_error(lk, method = "ajk_bc", vcov = V)

  expect_gt(max(le$le$le_mu), 0)
  expect_equal(le$le$lebc_mu, le$le$le_mu, tolerance = 1e-8)
  expect_equal(le$le$tebc_mu, le$le$le_mu, tolerance = 1e-8)
  expect_identical(le$le$method, rep("ajk_bc", nrow(le$le)))
})

test_that("analytic bias-correction derivatives match the numeric ones", {
  # Unbalanced presence patterns (one item missing in group 1, another in
  # group 3) and a non-diagonal item-parameter covariance make this
  # comparison sensitive to column-mapping errors in the analytic H_gamma.
  ip <- make_identified_ipars(c(0, 0.3, 0.6), c(1, 1.1, 1.2), I = 6)
  ip$b[ip$group == 2 & ip$item == "I01"] <- ip$b[ip$group == 2 & ip$item == "I01"] + 0.5
  ip$a[ip$group == 3 & ip$item == "I02"] <- ip$a[ip$group == 3 & ip$item == "I02"] * 1.2
  ip <- ip[!(ip$group == 1 & ip$item == "I06"), , drop = FALSE]
  ip <- ip[!(ip$group == 3 & ip$item == "I05"), , drop = FALSE]

  links <- list(link_chain(ip, method = "mgm"),
                link_chain(ip, method = "mm"),
                link_joint(ip, method = "phl"),
                link_joint(ip, method = "phl", item_weights = "inverse_admin"))
  for (lk in links) {
    p <- 2 * nrow(lk$ipars)
    set.seed(99)
    R <- matrix(stats::rnorm(p * p), p, p)
    V <- 0.0005 * (tcrossprod(R) / p + diag(p))
    an <- le_jackknife_bc(lk, kind = "jackknife", vcov = V,
                          derivatives = "analytic")
    nu <- le_jackknife_bc(lk, kind = "jackknife", vcov = V,
                          derivatives = "numeric")
    for (cc in c("le_mu", "le_sigma", "se_mu", "se_sigma",
                 "lebc_mu", "lebc_sigma", "te_mu", "te_sigma",
                 "tebc_mu", "tebc_sigma")) {
      expect_equal(an[[cc]], nu[[cc]], tolerance = 1e-6)
    }
  }
})

test_that("numeric bias-correction derivatives converge to the analytic ones", {
  # Central finite differences have error O(h^2), so shrinking h by a
  # factor of 10 must shrink the deviation from the analytic values by
  # roughly a factor of 100. This identifies the analytic derivatives as
  # the limit of the numeric approximation.
  ip <- make_identified_ipars(c(0, 0.3, 0.6), c(1, 1.1, 1.2), I = 6)
  ip$b[ip$group == 2 & ip$item == "I01"] <- ip$b[ip$group == 2 & ip$item == "I01"] + 0.5
  ip$a[ip$group == 3 & ip$item == "I02"] <- ip$a[ip$group == 3 & ip$item == "I02"] * 1.2

  cols <- c("se_mu", "se_sigma", "lebc_mu", "lebc_sigma",
            "tebc_mu", "tebc_sigma")
  for (lk in list(link_chain(ip, method = "mgm"),
                  link_joint(ip, method = "phl", item_weights = "inverse_admin"))) {
    V <- diag(0.001, nrow = 2 * nrow(lk$ipars))
    an <- le_jackknife_bc(lk, kind = "jackknife", vcov = V,
                          derivatives = "analytic")
    err <- vapply(c(1e-2, 1e-3, 1e-4), function(h) {
      nu <- le_jackknife_bc(lk, kind = "jackknife", vcov = V, h = h,
                            derivatives = "numeric")
      max(vapply(cols, function(cc) max(abs(an[[cc]] - nu[[cc]]),
                                        na.rm = TRUE), numeric(1)))
    }, numeric(1))
    expect_lt(err[2], err[1] / 20)
    expect_lt(err[3], err[2] / 20)
  }
})

test_that("bias-corrected AJK linking error equals the analytic AJK", {
  ip <- make_identified_ipars(c(0, 0.3, 0.6), c(1, 1.1, 1.2), I = 6)
  ip$b[ip$group == 2 & ip$item == "I01"] <- ip$b[ip$group == 2 & ip$item == "I01"] + 0.5
  ip$a[ip$group == 3 & ip$item == "I02"] <- ip$a[ip$group == 3 & ip$item == "I02"] * 1.2

  for (lk in list(link_chain(ip, method = "mgm"),
                  link_joint(ip, method = "phl", item_weights = "inverse_admin"))) {
    V <- diag(0.001, nrow = 2 * nrow(lk$ipars))
    bc <- linking_error(lk, method = "ajk_bc", vcov = V)
    ajk <- linking_error(lk, method = "ajk")

    expect_equal(bc$le$le_mu, ajk$le$le_mu, tolerance = 1e-10)
    expect_equal(bc$le$le_sigma, ajk$le$le_sigma, tolerance = 1e-10)
    expect_gt(max(bc$le$se_mu), 0)
    expect_true(all(bc$le$lebc_mu <= bc$le$le_mu + 1e-12))
    expect_true(all(bc$le$lebc_sigma <= bc$le$le_sigma + 1e-12))
  }
})

test_that("sandwich linking error is zero for identified drift-free chain links", {
  ip <- make_identified_ipars(c(0, 0.3, 0.6), c(1, 1.1, 1.2), I = 12)

  for (link_method in c("mgm", "mm")) {
    lk <- link_chain(ip, method = link_method)
    esw <- linking_error(lk, method = "sandwich_esw")
    osw <- linking_error(lk, method = "sandwich_osw")
    bosw <- linking_error(lk, method = "sandwich_bosw")

    expect_named(osw$le, c("group", "le_mu", "le_sigma", "method", "n_reps"))
    expect_lt(max(abs(esw$le$le_mu), na.rm = TRUE), 1e-8)
    expect_lt(max(abs(esw$le$le_sigma), na.rm = TRUE), 1e-8)
    expect_lt(max(abs(osw$le$le_mu), na.rm = TRUE), 1e-8)
    expect_lt(max(abs(osw$le$le_sigma), na.rm = TRUE), 1e-8)
    expect_lt(max(abs(bosw$le$le_mu), na.rm = TRUE), 1e-8)
    expect_lt(max(abs(bosw$le$le_sigma), na.rm = TRUE), 1e-8)
    expect_identical(esw$le$method, rep("sandwich_esw", nrow(esw$le)))
    expect_identical(osw$le$method, rep("sandwich_osw", nrow(osw$le)))
    expect_identical(bosw$le$method, rep("sandwich_bosw", nrow(bosw$le)))
  }
})

test_that("sandwich linking error is zero for identified drift-free PHL links", {
  ip <- make_identified_ipars(c(0, 0.3, 0.6), c(1, 1.1, 1.2), I = 12)

  for (iw in c("uniform", "inverse_admin")) {
    lk <- link_joint(ip, method = "phl", item_weights = iw)
    esw <- linking_error(lk, method = "sandwich_esw")
    le <- linking_error(lk, method = "sandwich_osw")

    expect_lt(max(abs(esw$le$le_mu), na.rm = TRUE), 1e-8)
    expect_lt(max(abs(esw$le$le_sigma), na.rm = TRUE), 1e-8)
    expect_lt(max(abs(le$le$le_mu), na.rm = TRUE), 1e-8)
    expect_lt(max(abs(le$le$le_sigma), na.rm = TRUE), 1e-8)
    expect_identical(esw$le$method, rep("sandwich_esw", nrow(esw$le)))
    expect_identical(le$le$method, rep("sandwich_osw", nrow(le$le)))
  }
})

test_that("sandwich linking error reacts to item-level drift", {
  ip <- make_identified_ipars(c(0, 0.3), c(1, 1.1), I = 12)
  ip$b[ip$group == 2 & ip$item == "I01"] <- ip$b[ip$group == 2 & ip$item == "I01"] + 0.5
  lk <- link_chain(ip, method = "mgm")
  esw <- linking_error(lk, method = "sandwich_esw")
  osw <- linking_error(lk, method = "sandwich_osw")
  bosw <- linking_error(lk, method = "sandwich_bosw")

  expect_gt(esw$le$le_mu[esw$le$group == 2], 0)
  expect_gt(osw$le$le_mu[osw$le$group == 2], 0)
  expect_gt(bosw$le$le_mu[bosw$le$group == 2],
            osw$le$le_mu[osw$le$group == 2])
  expect_true(all(osw$le$n_reps == 12))
})

test_that("expected sandwich reproduces Robitzsch two-group MGM formula", {
  I <- 12
  item <- sprintf("I%02d", seq_len(I))
  base_a <- rep(c(0.73, 1.25, 1.20, 1.47, 0.97, 1.38, 1.05, 1.14, 1.15, 0.67),
                length.out = I)
  base_b <- rep(c(-1.31, 1.44, -1.20, 0.10, 0.10, -0.74, 1.48, -0.61, 0.82,
                  -0.07), length.out = I)
  f <- seq(-0.20, 0.25, length.out = I)
  e <- c(0.40, -0.20, 0.10, 0.30, -0.10, 0.20,
         -0.30, 0.05, 0.15, -0.25, 0.35, -0.15)
  mu <- 0.3
  sigma <- 1.1
  ip <- rbind(
    data.frame(group = 1, item = item, a = base_a, b = base_b),
    data.frame(group = 2, item = item, a = base_a * sigma * exp(f),
               b = (base_b + e - mu) / sigma)
  )
  lk <- link_chain(ip, method = "mgm")
  esw <- linking_error(lk, method = "sandwich_esw")

  v1 <- log(ip$a[ip$group == 2]) - log(ip$a[ip$group == 1])
  sigma_est <- exp(mean(v1))
  v2 <- sigma_est * ip$b[ip$group == 2] - ip$b[ip$group == 1]
  mu_est <- -mean(v2)
  tau_a <- stats::sd(v1)
  tau_b <- stats::sd(v2)
  tau_ab <- stats::cov(v1, v2)
  h1 <- mu_est - mean(ip$b[ip$group == 1])
  expected_mu <- sqrt(tau_b^2 / I + h1^2 * tau_a^2 / I -
                        2 * h1 * tau_ab / I)
  expected_sigma <- sqrt(sigma_est^2 * tau_a^2 / I)

  expect_equal(esw$le$le_mu[esw$le$group == 2], expected_mu,
               tolerance = 1e-10)
  expect_equal(esw$le$le_sigma[esw$le$group == 2], expected_sigma,
               tolerance = 1e-10)
})

test_that("analytic observed sandwich matches observed sandwich for supported links", {
  ip <- make_identified_ipars(c(0, 0.3, 0.6), c(1, 1.1, 1.2), I = 12)
  ip$b[ip$group == 2 & ip$item == "I01"] <- ip$b[ip$group == 2 & ip$item == "I01"] + 0.5
  ip$a[ip$group == 2 & ip$item == "I02"] <- ip$a[ip$group == 2 & ip$item == "I02"] * 1.2
  ip$b[ip$group == 3 & ip$item == "I03"] <- ip$b[ip$group == 3 & ip$item == "I03"] - 0.4
  ip$a[ip$group == 3 & ip$item == "I04"] <- ip$a[ip$group == 3 & ip$item == "I04"] * 0.85

  for (link_method in c("mgm", "mm")) {
    lk <- link_chain(ip, method = link_method)
    osw <- linking_error(lk, method = "sandwich_osw")
    an <- linking_error(lk, method = "sandwich_osw_analytic")

    expect_equal(an$le$le_mu, osw$le$le_mu, tolerance = 1e-8)
    expect_equal(an$le$le_sigma, osw$le$le_sigma, tolerance = 1e-8)
    expect_identical(an$le$method, rep("sandwich_osw_analytic", nrow(an$le)))
  }

  for (iw in c("uniform", "inverse_admin")) {
    lk <- link_joint(ip, method = "phl", item_weights = iw)
    osw <- linking_error(lk, method = "sandwich_osw")
    an <- linking_error(lk, method = "sandwich_osw_analytic")

    expect_equal(an$le$le_mu, osw$le$le_mu, tolerance = 1e-8)
    expect_equal(an$le$le_sigma, osw$le$le_sigma, tolerance = 1e-8)
    expect_identical(an$le$method, rep("sandwich_osw_analytic", nrow(an$le)))
  }
})

test_that("linking_error validates inputs and reserved methods", {
  expect_error(linking_error(list()), "irtlink_link")
  ip <- make_identified_ipars(c(0, 0.3), c(1, 1.1), I = 3)
  lk <- link_chain(ip, method = "mgm", anchor = "I01")
  expect_error(linking_error(lk), "at least two")
  unsupported_chain <- link_chain(ip, method = "mgm")
  unsupported_chain$method <- "haebara"
  unsupported_restricted <- link_joint(ip, method = "phl")
  unsupported_restricted$approach <- "joint_restricted"
  expect_error(linking_error(unsupported_chain, method = "sandwich_osw"),
               "chain links")
  expect_error(linking_error(unsupported_restricted, method = "sandwich_osw"),
               "approach")
  expect_error(linking_error(unsupported_chain, method = "sandwich_esw"),
               "chain links")
  expect_error(linking_error(unsupported_restricted, method = "sandwich_esw"),
               "approach")
  # chain Haebara and SL now have an AJK; Haberman remains rejected
  unsupported_chain_hab <- link_chain(ip, method = "mgm")
  unsupported_chain_hab$method <- "haberman"
  expect_error(linking_error(unsupported_chain_hab, method = "ajk"),
               "chain links")
  expect_error(linking_error(unsupported_chain, method = "jackknife_full"),
               "should be one of")
  expect_error(linking_error(unsupported_chain, method = "ajk_full"),
               "should be one of")
  expect_error(linking_error(link_chain(ip, method = "mgm"),
                             method = "ajk_bc"),
               "needs `vcov`")
  expect_error(linking_error(unsupported_restricted, method = "ajk"),
               "approach")
  expect_error(linking_error(unsupported_restricted, method = "ajk_bc",
                             vcov = diag(2 * nrow(unsupported_restricted$ipars))),
               "approach")
  expect_error(linking_error(unsupported_chain, method = "sandwich_osw_analytic"),
               "chain links")
  expect_error(linking_error(unsupported_restricted,
                             method = "sandwich_osw_analytic"),
               "approach")
})
