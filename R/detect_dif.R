# R/detect_dif.R

detect_dif <- function(calib, method = c("rmsd", "lrt"),
                       rmsd_method = c("parameter", "data", "both"),
                       cutoff_type = c("fixed", "data_driven"),
                       cutoff = 0.05, tau = 2.7,
                       flag_on = c("parameter", "data"),
                       link_method = "mgm",
                       procedure = c("one_step", "iterative_forward"),
                       max_iter = 7, min_anchor = 3,
                       alpha = 0.05,
                       alpha_adjust = c("none", "bonferroni"),
                       dif_type = c("both", "uniform", "nonuniform"),
                       ...) {
  if (!inherits(calib, "irtlink_calib")) {
    stop("`calib` must be an `irtlink_calib` object (from calibrate()).",
         call. = FALSE)
  }
  if (calib$model %in% c("GPCM", "PCM")) {
    stop("DIF detection for polytomous models is not yet available.",
         call. = FALSE)
  }
  if (!is.null(calib$pweights)) {
    warning("DIF detection refits the groups without the person ",
            "weights of the calibration.", call. = FALSE)
  }
  flag_on_was_passed <- !missing(flag_on)
  dif_type_passed <- !missing(dif_type)
  method <- match.arg(method)
  rmsd_method <- match.arg(rmsd_method)
  cutoff_type <- match.arg(cutoff_type)
  flag_on <- match.arg(flag_on)
  procedure <- match.arg(procedure)
  alpha_adjust <- match.arg(alpha_adjust)
  dif_type <- match.arg(dif_type)
  if (rmsd_method != "both" && flag_on_was_passed) {
    warning("`flag_on` is ignored when rmsd_method is not \"both\".",
            call. = FALSE)
  }
  if (cutoff_type == "data_driven" && length(tau) != 1L) {
    stop("`tau` must be a single value (cutoff sensitivity vectors are ",
         "only supported for fixed cutoffs).", call. = FALSE)
  }
  if (method == "lrt") {
    if (procedure == "iterative_forward")
      warning("LRT is one-step only; `procedure` is ignored.", call. = FALSE)
    res <- detect_lrt(calib, alpha = alpha, alpha_adjust = alpha_adjust,
                      dif_type = dif_type, min_anchor = min_anchor, ...)
    res$method <- method
    res$rmsd_method <- NA_character_
    res$cutoff_type <- NA_character_
    res$cutoff <- NA_real_
    res$tau <- NA_real_
    res$alpha <- alpha
    res$alpha_adjust <- alpha_adjust
    res$dif_type <- dif_type
    res$min_anchor <- min_anchor
    res$procedure <- "one_step"
    res$call <- match.call()
    return(structure(res, class = "irtlink_dif"))
  }
  if (dif_type_passed && dif_type != "both")
    warning("dif_type applies only to method = \"lrt\"; ignored for RMSD.",
            call. = FALSE)
  driver <- if (rmsd_method == "both") flag_on else rmsd_method
  if (rmsd_method == "both" && procedure == "iterative_forward") {
    warning("rmsd_method = \"both\" with iterative_forward runs only the ",
            "flag_on estimator's iterative purification; ",
            "the other estimator's table is NULL.",
            call. = FALSE)
  }

  if (procedure == "one_step") {
    res <- detect_dif_onestep(calib, rmsd_method, driver, cutoff_type,
                              cutoff, tau, link_method, min_anchor, ...)
  } else if (driver == "data") {
    res <- detect_dif_iterative_data(calib, cutoff_type, cutoff, tau,
                                     max_iter, min_anchor, ...)
  } else {
    res <- detect_dif_iterative(calib, cutoff_type, cutoff, tau,
                                link_method, max_iter, min_anchor, ...)
  }
  res$method <- method
  res$rmsd_method <- rmsd_method
  res$cutoff_type <- cutoff_type
  res$cutoff <- cutoff
  res$tau <- tau
  res$min_anchor <- min_anchor
  res$procedure <- procedure
  res$lrt <- NULL
  res$alpha <- NA_real_
  res$alpha_adjust <- NA_character_
  res$dif_type <- NA_character_
  res$call <- match.call()
  structure(res, class = "irtlink_dif")
}

# requested RMSD tables, driver points at the flagging table
compute_rmsd_tables <- function(calib, rmsd_method, driver,
                                link_method, anchor, ...) {
  tab_p <- if (rmsd_method %in% c("parameter", "both")) {
    rmsd_parameter(calib$ipars, link_method = link_method, anchor = anchor)
  } else NULL
  tab_d <- if (rmsd_method %in% c("data", "both")) {
    rmsd_data(calib, ...)
  } else NULL
  list(parameter = tab_p, data = tab_d,
       driver = if (driver == "data") tab_d else tab_p)
}

detect_dif_onestep <- function(calib, rmsd_method, driver, cutoff_type,
                               cutoff, tau, link_method, min_anchor, ...) {
  tabs <- compute_rmsd_tables(calib, rmsd_method, driver,
                              link_method, anchor = NULL, ...)
  tab <- tabs$driver
  linkable <- tab$item
  cuts <- if (cutoff_type == "fixed") cutoff else tau
  by_cut <- lapply(cuts, function(cv) {
    if (cutoff_type == "fixed") {
      dif_flag_fixed(tab, cutoff = cv, min_anchor = min_anchor)
    } else {
      wcols <- grep("^g[0-9]+$", names(tab), value = TRUE)
      dif_flag_mad(tab, tau = cv, min_anchor = min_anchor, group_cols = wcols)
    }
  })
  names(by_cut) <- as.character(cuts)
  flagged <- by_cut[[1]]
  list(
    rmsd = tab, rmsd_parameter = tabs$parameter, rmsd_data = tabs$data,
    flagged = flagged, anchor = setdiff(linkable, flagged),
    flagged_by_cutoff = by_cut, linkable = linkable, history = NULL,
    converged = TRUE
  )
}

# forward purification on the parameter-based RMSD, re-link each
# round, never retest flagged items
detect_dif_iterative <- function(calib, cutoff_type, cutoff, tau,
                                 link_method, max_iter, min_anchor, ...) {
  ip <- as_ipars(calib$ipars)
  presence <- table(ip$item)
  linkable <- names(presence)[presence >= 2]
  anchor <- linkable
  flagged <- character(0)
  history <- list()
  cut1 <- if (cutoff_type == "fixed") cutoff[1] else tau

  if (cutoff_type == "fixed" && length(cutoff) > 1L) {
    warning("cutoff vectors are only supported for procedure = ",
            "\"one_step\"; using cutoff[1] = ", cutoff[1],
            " for iterative_forward.", call. = FALSE)
  }
  if (length(linkable) == 0L) {
    empty_tab <- data.frame(item = character(0), max = numeric(0))
    return(list(
      rmsd = empty_tab, rmsd_parameter = empty_tab, rmsd_data = NULL,
      flagged = character(0), anchor = character(0),
      flagged_by_cutoff = stats::setNames(list(character(0)),
                                          as.character(cut1)),
      linkable = character(0), history = list(),
      converged = TRUE
    ))
  }

  hit_min_anchor <- FALSE
  for (iter in seq_len(max_iter)) {
    tab <- rmsd_parameter(ip, link_method = link_method, anchor = anchor)
    cand <- tab[tab$item %in% anchor, , drop = FALSE]
    if (cutoff_type == "fixed") {
      new_flags <- dif_flag_fixed(cand, cutoff = cut1, min_anchor = min_anchor)
    } else {
      wcols <- grep("^g[0-9]+$", names(cand), value = TRUE)
      new_flags <- dif_flag_mad(cand, tau = cut1, min_anchor = min_anchor,
                                group_cols = wcols)
    }
    history[[iter]] <- list(anchor = anchor, rmsd = tab, new_flags = new_flags)
    if (length(new_flags) == 0) break
    flagged <- c(flagged, new_flags)
    anchor <- setdiff(anchor, new_flags)
    if (length(anchor) <= min_anchor) {
      hit_min_anchor <- TRUE
      break
    }
  }

  last_new <- history[[length(history)]]$new_flags
  converged <- length(last_new) == 0L && !hit_min_anchor
  if (hit_min_anchor)
    warning("iterative_forward: stopped early because the anchor would drop ",
            "to min_anchor (", min_anchor, "); inspect $history.",
            call. = FALSE)
  if (!converged && !hit_min_anchor && length(history) == max_iter) {
    warning("iterative_forward: max_iter (", max_iter, ") reached before ",
            "convergence; purification may be incomplete. Increase max_iter ",
            "or inspect $history.", call. = FALSE)
  }

  list(
    rmsd = history[[length(history)]]$rmsd,
    rmsd_parameter = history[[length(history)]]$rmsd, rmsd_data = NULL,
    flagged = flagged, anchor = anchor,
    flagged_by_cutoff = stats::setNames(list(flagged), as.character(cut1)),
    linkable = linkable, history = history,
    converged = converged
  )
}
