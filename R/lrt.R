# R/lrt.R
# likelihood ratio test for DIF/IPD, for each tested item one model
# that frees this item across groups

lrt_df <- function(dif_type, n_groups) {
  base <- n_groups - 1L
  if (dif_type == "both") 2L * base else base
}

# per-item LRT statistics over the common items
lrt_item_deviances <- function(calib, dif_type = c("both", "uniform",
                                                   "nonuniform"),
                               verbose = FALSE, control = list(), ...) {
  dif_type <- match.arg(dif_type)
  rpw <- responses_per_group(calib)            # enforces separate
  # deviance differences need a tighter EM convergence than
  # parameter estimation
  ctrl <- utils::modifyList(list(conv = 1e-7, maxit = 2000L), control)
  ctrl$verbose <- verbose
  base_dev <- em_fit_multi(rpw, model = calib$model,
                           control = ctrl)$deviance

  presence <- table(unlist(lapply(rpw, colnames)))
  common <- names(presence)[presence >= 2]

  out <- lapply(common, function(it) {
    # df from the number of groups the item appears in
    n_present <- as.integer(presence[[it]])
    df_it <- lrt_df(dif_type, n_present)
    aug <- em_fit_multi(rpw, model = calib$model, control = ctrl,
                        free_items = if (dif_type %in%
                                         c("both", "uniform")) it
                                     else NULL,
                        free_slope_items = if (dif_type %in%
                                               c("both", "nonuniform")) it
                                           else NULL)$deviance
    gap <- base_dev - aug
    # a clearly negative gap signals a convergence problem
    if (gap < -1e-2)
      warning("LRT: augmented model fits worse than the base for item ", it,
              " (deviance gap ", round(gap, 3), "); chi-square set to 0 ",
              "(check convergence).", call. = FALSE)
    chisq <- max(0, gap)
    data.frame(item = it, chisq = chisq, df = df_it,
               p = stats::pchisq(chisq, df_it, lower.tail = FALSE),
               stringsAsFactors = FALSE)
  })
  do.call(rbind, out)
}

# flag items from the per-item LRT, optional Bonferroni adjustment,
# at least min_anchor items remain
detect_lrt <- function(calib, alpha, alpha_adjust, dif_type, min_anchor, ...) {
  lt <- lrt_item_deviances(calib, dif_type = dif_type, ...)
  p <- lt$p
  if (alpha_adjust == "bonferroni") p <- stats::p.adjust(p, method = "bonferroni")
  lt$p_adj <- p
  above <- lt$item[!is.na(p) & p < alpha]
  stat <- stats::setNames(lt$chisq, lt$item)   # tie-break by chisq
  flagged <- dif_tiebreak(stat, above, min_anchor = min_anchor,
                          total = nrow(lt))
  list(
    rmsd = NULL, rmsd_parameter = NULL, rmsd_data = NULL, lrt = lt,
    flagged = flagged, anchor = setdiff(lt$item, flagged),
    flagged_by_cutoff = stats::setNames(list(flagged), as.character(alpha)),
    linkable = lt$item, history = NULL, converged = TRUE
  )
}
