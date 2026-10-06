# R/estimate_trend.R

estimate_trend <- function(data_list,
                           calibration = c("separate", "concurrent",
                                           "fixed", "regularized"),
                           approach = c("chain", "joint", "joint_restricted"),
                           link = c("haberman", "haebara", "sl", "mm", "mgm"),
                           variant = c("simultaneous", "pairwise"),
                           pow = 2, weights = "uniform", use_intercepts = TRUE,
                           item_weights = c("inverse_admin", "uniform"),
                           dif = c("none", "purify"),
                           le = c("none", "jackknife", "jackknife_bc",
                                   "ajk", "ajk_bc",
                                   "sandwich_esw", "sandwich_osw",
                                   "sandwich_bosw", "sandwich_esw_bc",
                                   "sandwich_osw_bc", "sandwich_bosw_bc",
                                   "sandwich_osw_analytic"),
                           le_vcov = NULL, le_h = 1e-5, ref = NULL,
                           model = c("2PL", "1PL", "GPCM", "PCM"),
                           eps = 0.001, ...) {
  calibration <- match.arg(calibration)
  approach <- match.arg(approach)
  link <- match.arg(link)
  variant <- match.arg(variant)
  item_weights <- match.arg(item_weights)
  dif <- match.arg(dif)
  le <- match.arg(le)
  model <- match.arg(model)
  if (length(data_list) < 2L)
    stop("estimate_trend needs at least two groups.", call. = FALSE)
  if (model %in% c("GPCM", "PCM")) {
    if (calibration != "separate")
      stop("model = \"", model, "\" currently supports only ",
           "calibration = \"separate\".", call. = FALSE)
    if (dif != "none")
      stop("DIF purification for polytomous models is not yet ",
           "available.", call. = FALSE)
  }
  if (variant == "pairwise") {
    if (approach == "chain")
      stop("variant = \"pairwise\" needs approach = \"joint\" or ",
           "\"joint_restricted\".", call. = FALSE)
    if (link %in% c("mm", "mgm"))
      stop("link = \"", link, "\" has no pairwise variant.", call. = FALSE)
    link <- switch(link, haberman = "phl", haebara = "haebara_pw",
                   sl = "sl_pw")
  }
  validate_trend_matrix(calibration, approach, link)

  le_bc_methods <- c("jackknife_bc", "ajk_bc", "sandwich_esw_bc",
                     "sandwich_osw_bc", "sandwich_bosw_bc")
  auto_vcov <- is.null(le_vcov) && le %in% le_bc_methods

  res <- switch(calibration,
    separate = trend_separate(data_list, approach, link, pow, weights,
                              use_intercepts, item_weights, dif, model,
                              keep_vcov = auto_vcov, ref = ref, ...),
    concurrent = trend_concurrent(data_list, approach, dif, model, ...),
    fixed = trend_fixed(data_list, approach, dif, model, ...),
    regularized = trend_regularized(data_list, approach, dif, model, eps, ...))
  if (calibration != "separate") {
    res$trend <- rebase_trend(res$trend, ref)
  }

  if (le != "none") {
    if (is.null(res$link_obj))
      stop("le = \"", le, "\" is currently available only for ",
           "separate-calibration linking results.", call. = FALSE)
    vcov_arg <- le_vcov
    if (is.null(vcov_arg) && le %in% le_bc_methods) {
      vcov_arg <- res$calib$vcov_ipars
      if (is.null(vcov_arg)) {
        stop("le = \"", le, "\" needs item-parameter covariance matrices ",
             "from the calibration, or pass `le_vcov` explicitly.",
             call. = FALSE)
      }
      message("le = \"", le, "\" uses the per-group covariance matrices ",
              "from the calibration, which treat the person samples of ",
              "the groups as independent. Pass `le_vcov` (for example ",
              "from dependent_vcov()) when the samples are dependent.")
    }
    res$link_obj <- linking_error(res$link_obj, method = le, vcov = vcov_arg,
                                  h = le_h)
    res$le_obj <- res$link_obj$le
  } else {
    res$le_obj <- NULL
  }

  res$trend <- finalize_trend_table(res$trend, res$le_obj)
  res$calibration <- calibration
  res$approach <- approach
  res$link <- if (calibration == "separate") link else NA_character_
  res$call <- match.call()
  structure(res, class = "irtlink_trend")
}

# fixed (FIPC) is sequential and ignores `approach`
validate_trend_matrix <- function(calibration, approach, link) {
  if (calibration == "regularized" && approach != "chain")
    stop("regularized (SBIC) linking with approach = \"", approach,
         "\" is not yet supported; use approach = \"chain\".", call. = FALSE)
  if (calibration == "separate" && link %in% c("mm", "mgm") &&
      approach != "chain")
    stop("link = \"", link, "\" is available for approach = \"chain\" ",
         "only. Use approach = \"chain\" or link = \"haberman\", ",
         "\"haebara\", or \"sl\".", call. = FALSE)
}

# add the error placeholder columns to the trend table
finalize_trend_table <- function(trend, le = NULL) {
  trend$se <- NA_real_
  trend$le <- NA_real_
  trend$te <- NA_real_
  if (!is.null(le)) {
    idx <- match(trend$group, le$group)
    trend$le <- le$le_mu[idx]
    if ("se_mu" %in% names(le)) {
      trend$se <- le$se_mu[idx]
    }
    if ("te_mu" %in% names(le)) {
      trend$te <- le$te_mu[idx]
    }
    if ("lebc_mu" %in% names(le)) {
      trend$le_bc <- le$lebc_mu[idx]
    }
    if ("tebc_mu" %in% names(le)) {
      trend$te_bc <- le$tebc_mu[idx]
    }
  }
  trend
}

# split dots into detect_dif and calibrate argument lists
split_trend_dots <- function(dots) {
  detect_nms <- c("method", "rmsd_method", "cutoff_type", "cutoff", "tau",
                  "flag_on", "link_method", "procedure", "max_iter",
                  "min_anchor", "alpha", "alpha_adjust", "dif_type")
  list(detect = dots[intersect(names(dots), detect_nms)],
       calib  = dots[setdiff(names(dots), detect_nms)])
}

# --- separate path ----------------------------------------------------------
trend_separate <- function(data_list, approach, link, pow, weights,
                           use_intercepts, item_weights, dif, model,
                           keep_vcov = FALSE, ref = NULL, ...) {
  parts <- split_trend_dots(list(...))
  if (keep_vcov) {
    parts$calib$keep_vcov <- TRUE
  }
  cal <- do.call(calibrate, c(list(data_list, model = model,
                                   calibration = "separate"), parts$calib))
  dif_obj <- NULL
  anchor <- NULL
  if (dif == "purify") {
    dif_obj <- do.call(detect_dif, c(list(cal), parts$detect))
    anchor <- dif_obj$anchor
  }
  link_obj <- if (approach == "chain") {
    link_chain(cal, method = link, pow = pow, weights = weights,
               use_intercepts = use_intercepts, anchor = anchor, ref = ref)
  } else {
    link_joint(cal, method = link, pow = pow, use_intercepts = use_intercepts,
               weights = weights, item_weights = item_weights,
               restricted = (approach == "joint_restricted"), anchor = anchor,
               ref = ref)
  }
  list(trend = link_obj$trend, calib = cal, dif = dif_obj, link_obj = link_obj)
}

# --- fixed / regularized paths (trend comes from calibration) ---------------
# fixed (FIPC) is sequential and ignores `approach`, regularized
# is chain-only
trend_fixed <- function(data_list, approach, dif, model, ...) {
  if (dif == "purify")
    stop("dif = \"purify\" is not supported for fixed (FIPC) calibration ",
         "(no partial-invariance mechanism). Use calibration = ",
         "\"concurrent\"/\"separate\", or dif = \"none\".", call. = FALSE)
  parts <- split_trend_dots(list(...))
  cal <- do.call(calibrate, c(list(data_list, model = model,
                                   calibration = "fixed"),
                              parts$calib))
  list(trend = cal$trend, calib = cal, dif = NULL, link_obj = NULL)
}

trend_regularized <- function(data_list, approach, dif, model, eps, ...) {
  if (dif == "purify")
    stop("dif = \"purify\" is redundant for regularized (SBIC) calibration: ",
         "SBIC already regularizes differential item functioning. Use dif = \"none\".",
         call. = FALSE)
  parts <- split_trend_dots(list(...))
  cal <- do.call(calibrate, c(list(data_list, model = model,
                                   calibration = "regularized",
                                   eps = eps), parts$calib))
  list(trend = cal$trend, calib = cal, dif = NULL, link_obj = NULL)
}

# --- concurrent path --------------------------------------------------------
trend_concurrent <- function(data_list, approach, dif, model, ...) {
  parts <- split_trend_dots(list(...))
  # exact [[, `$` would partial-match another calib argument to `engine`
  engine <- parts$calib[["engine"]]
  if (is.null(engine)) engine <- "em"
  engine <- match.arg(engine, "em")
  parts$calib[["engine"]] <- NULL
  if (dif == "purify" && approach != "joint") {
    stop("dif = \"purify\" with concurrent calibration is only supported ",
         "for approach = \"joint\" in this version.", call. = FALSE)
  }
  free_items <- NULL
  dif_obj <- NULL
  if (dif == "purify") {
    sep <- do.call(calibrate, c(list(data_list, model = model,
                                     calibration = "separate",
                                     engine = engine), parts$calib))
    dif_obj <- do.call(detect_dif, c(list(sep), parts$detect))
    free_items <- dif_obj$flagged
    # nothing flagged, fall back to the full-invariance concurrent fit
    if (length(free_items) == 0) free_items <- NULL
  }
  fit <- switch(approach,
    joint = do.call(concurrent_window,
                    c(list(data_list, model = model,
                           free_items = free_items, engine = engine),
                      parts$calib)),
    chain = do.call(fit_concurrent_chain,
                    c(list(data_list, model = model, engine = engine),
                      parts$calib)),
    joint_restricted = do.call(fit_concurrent_restricted,
                               c(list(data_list, model = model,
                                      engine = engine), parts$calib)))
  cal <- structure(
    list(ipars = fit$ipars, model = model, calibration = "concurrent",
         engine = engine, n_groups = length(data_list),
         N = vapply(data_list, nrow, integer(1)),
         trend = fit$trend, converged = fit$converged),
    class = "irtlink_calib")
  list(trend = fit$trend, calib = cal, dif = dif_obj, link_obj = NULL)
}
