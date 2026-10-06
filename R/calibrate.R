# R/calibrate.R

calibrate <- function(data_list, model = c("2PL", "1PL", "3PL",
                                           "GPCM", "PCM"),
                      calibration = c("separate", "concurrent",
                                      "fixed", "regularized"),
                      engine = "em",
                      free_items = NULL,
                      eps = 0.001,
                      pweights = NULL,
                      keep_models = TRUE, keep_vcov = FALSE,
                      verbose = FALSE, ...) {
  model <- match.arg(model)
  engine <- match.arg(engine)
  calibration <- match.arg(calibration)
  if (!is.null(free_items) && calibration != "concurrent")
    stop("free_items requires calibration = \"concurrent\".", call. = FALSE)
  if (model %in% c("3PL", "GPCM", "PCM") && calibration != "separate")
    stop("model = \"", model, "\" currently supports only calibration = ",
         "\"separate\".", call. = FALSE)
  design <- linking_design(data_list)
  if (model %in% c("GPCM", "PCM")) {
    check_polytomous(data_list)
  } else {
    check_dichotomous(data_list)
  }
  check_pweights(pweights, data_list)
  if (any(design$successive_common == 0)) {
    stop("Successive groups share no common items; linking is not possible.",
         call. = FALSE)
  }
  if (any(design$successive_common < 5)) {
    warning("Fewer than 5 common items between some successive groups; ",
            "linking may be unstable.", call. = FALSE)
  }

  trend    <- NULL
  eps_used <- NULL
  vcov_ipars <- NULL
  item_parameters <- NULL
  g_estimates <- NULL

  if (calibration == "separate") {
    sep <- lapply(seq_along(data_list), function(t)
      fit_em_separate(data_list[[t]], model = model, group = t,
                      verbose = verbose,
                      pweights = if (is.null(pweights)) NULL
                                 else pweights[[t]], ...))
    fits <- lapply(sep, `[[`, "model")
    ipars <- do.call(rbind, lapply(sep, `[[`, "ipars"))
    item_parameters <- add_intercept_columns(ipars, "difficulty")
    converged <- vapply(sep, `[[`, logical(1), "converged")
    if (keep_vcov) {
      vcov_ipars <- lapply(fits, em_vcov_ipars)
      names(vcov_ipars) <- as.character(seq_along(vcov_ipars))
    }
  } else if (calibration == "concurrent") {
    fc <- fit_em_concurrent(data_list, model = model,
                            free_items = free_items,
                            verbose = verbose, pweights = pweights, ...)
    ipars <- fc$ipars
    item_parameters <- add_intercept_columns(ipars, "difficulty")
    fits <- list(fc$model)
    converged <- stats::setNames(fc$converged, "concurrent_fit")
    trend <- fc$trend
    g_estimates <- fc$g
  } else if (calibration == "fixed") {
    fc <- fit_em_fipc(data_list, model = model, verbose = verbose,
                      pweights = pweights, ...)
    ipars <- fc$ipars
    item_parameters <- add_intercept_columns(ipars, "difficulty")
    fits <- NULL          # FIPC is sequential; no single representative fit
    converged <- stats::setNames(fc$converged, "fipc_fit")
    trend <- fc$trend
  } else if (calibration == "regularized") {
    # one fit for a single eps, a warm-started grid for a vector
    fc <- fit_em_sbic(data_list, model = model,
                      eps = sort(unique(eps), decreasing = TRUE),
                      verbose = verbose, pweights = pweights, ...)
    ipars     <- fc$ipars
    item_parameters <- add_intercept_columns(ipars, "difficulty")
    fits      <- NULL
    converged <- stats::setNames(fc$converged, "sbic_fit")
    trend     <- fc$trend
    eps_used  <- fc$eps
    g_estimates <- fc$g
  } else {
    stop("calibration = \"", calibration, "\" is not yet implemented.",
         call. = FALSE)
  }
  if (keep_vcov && is.null(vcov_ipars)) {
    warning("keep_vcov = TRUE is currently supported only for ",
            "calibration = \"separate\"; vcov_ipars is NULL.",
            call. = FALSE)
  }

  rownames(ipars) <- NULL

  # normalize to the canonical ipars format
  ipars <- as_ipars(ipars)
  item_parameters <- normalize_item_parameters(item_parameters, ipars)

  structure(
    list(
      ipars = ipars,
      item_parameters = item_parameters,
      model = model,
      calibration = calibration,
      engine = engine,
      item_parameterization = "difficulty",
      n_groups = length(data_list),
      N = vapply(data_list, nrow, integer(1)),
      converged = converged,
      models = if (keep_models) fits else NULL,
      vcov_ipars = vcov_ipars,
      design = design,
      trend = trend,
      free_items = free_items,
      g_estimates = g_estimates,
      pweights = pweights,
      eps = eps_used,
      call = match.call()
    ),
    class = "irtlink_calib"
  )
}

add_intercept_columns <- function(ipars, estimated_parameterization = "difficulty") {
  out <- ipars
  out$d <- -out$a * out$b
  out$nu <- out$a * out$b
  out$estimated_parameterization <- estimated_parameterization
  out
}

normalize_item_parameters <- function(item_parameters, ipars) {
  if (is.null(item_parameters)) {
    return(add_intercept_columns(ipars, "difficulty"))
  }
  required <- c("group", "item", "a", "b", "d", "nu",
                "estimated_parameterization")
  missing_cols <- setdiff(required, names(item_parameters))
  if (length(missing_cols) > 0) {
    stop("Internal error: item_parameters is missing columns: ",
         paste(missing_cols, collapse = ", "), call. = FALSE)
  }
  # keep the guessing and threshold columns of the 3PL and GPCM
  extras <- intersect(c("c", grep("^tau[0-9]+$", names(item_parameters),
                                  value = TRUE)),
                      names(item_parameters))
  item_parameters <- item_parameters[c(required, extras)]
  item_parameters$group <- as.integer(item_parameters$group)
  item_parameters$item <- as.character(item_parameters$item)
  item_parameters <- item_parameters[order(item_parameters$group,
                                           item_parameters$item), ]
  rownames(item_parameters) <- NULL
  item_parameters
}

check_dichotomous <- function(data_list) {
  for (t in seq_along(data_list)) {
    vals <- unlist(data_list[[t]], use.names = FALSE)
    vals <- vals[!is.na(vals)]
    if (!all(vals %in% c(0, 1))) {
      stop("Group ", t, " contains non-dichotomous responses; ",
           "only 0/1/NA are supported.", call. = FALSE)
    }
  }
  invisible(TRUE)
}

# person sampling weights, one positive vector per group
check_pweights <- function(pweights, data_list) {
  if (is.null(pweights)) return(invisible(TRUE))
  if (!is.list(pweights) || length(pweights) != length(data_list)) {
    stop("`pweights` must be a list with one numeric vector per group.",
         call. = FALSE)
  }
  for (t in seq_along(pweights)) {
    wt <- pweights[[t]]
    nt <- nrow(as.data.frame(data_list[[t]]))
    if (!is.numeric(wt) || length(wt) != nt || anyNA(wt) ||
        any(!is.finite(wt)) || any(wt <= 0)) {
      stop("`pweights[[", t, "]]` must hold ", nt,
           " positive finite values.", call. = FALSE)
    }
  }
  invisible(TRUE)
}

# integer categories 0..8 without gaps, every category observed
check_polytomous <- function(data_list) {
  for (t in seq_along(data_list)) {
    X <- as.matrix(as.data.frame(data_list[[t]]))
    vals <- X[!is.na(X)]
    if (!all(vals == floor(vals)) || any(vals < 0) || any(vals > 8)) {
      stop("Group ", t, " contains invalid responses; polytomous ",
           "models expect integer categories 0..8.", call. = FALSE)
    }
    for (j in seq_len(ncol(X))) {
      v <- X[, j]
      v <- v[!is.na(v)]
      if (!length(v)) next
      mx <- max(v)
      if (mx == 0) {
        stop("Item ", colnames(X)[j], " in group ", t,
             " has a single observed category.", call. = FALSE)
      }
      if (!all(seq(0, mx) %in% v)) {
        stop("Item ", colnames(X)[j], " in group ", t,
             " skips response categories; every category between 0 ",
             "and the item maximum must be observed.", call. = FALSE)
      }
    }
  }
  invisible(TRUE)
}
