# R/purify-data.R
# iterative item purification on the data-based RMSD, flagged items
# are split and refitted until no new flags appear
detect_dif_iterative_data <- function(calib, cutoff_type, cutoff, tau,
                                      max_iter, min_anchor, ...) {
  rpw <- responses_per_group(calib)           # enforces separate
  cut1 <- if (cutoff_type == "fixed") cutoff[1] else tau
  if (cutoff_type == "fixed" && length(cutoff) > 1L) {
    warning("cutoff vectors are only supported for procedure = ",
            "\"one_step\"; using cutoff[1] = ", cutoff[1],
            " for iterative_forward.", call. = FALSE)
  }

  flag_tab <- function(tab) {
    if (cutoff_type == "fixed") {
      dif_flag_fixed(tab, cutoff = cut1, min_anchor = min_anchor)
    } else {
      wcols <- grep("^g[0-9]+$", names(tab), value = TRUE)
      dif_flag_mad(tab, tau = cut1, min_anchor = min_anchor, group_cols = wcols)
    }
  }

  tab0 <- em_rmsd_data(calib)
  linkable <- tab0$item
  flagged <- flag_tab(tab0)
  # history rows carry anchor and split per round
  history <- list(list(anchor = linkable, split = character(0),
                       rmsd = tab0, new_flags = flagged))
  hit_min_anchor <- FALSE

  for (iter in seq_len(max_iter)) {
    if (length(flagged) == 0L) break
    if (length(setdiff(linkable, flagged)) <= min_anchor) {
      hit_min_anchor <- TRUE
      break
    }
    tab <- em_rmsd_data_split(rpw, split_items = flagged,
                              model = calib$model)
    cand <- tab[tab$item %in% setdiff(linkable, flagged), , drop = FALSE]
    new_flags <- flag_tab(cand)
    history[[length(history) + 1L]] <-
      list(anchor = setdiff(linkable, flagged), split = flagged,
           rmsd = tab, new_flags = new_flags)
    if (length(new_flags) == 0L) break
    flagged <- union(flagged, new_flags)
  }

  last_new <- history[[length(history)]]$new_flags
  converged <- length(last_new) == 0L && !hit_min_anchor
  if (!converged && hit_min_anchor) {
    warning("iterative_forward (data): stopped early because the anchor ",
            "would drop to min_anchor (", min_anchor, "); inspect $history.",
            call. = FALSE)
  }
  if (!converged && length(history) - 1L >= max_iter) {
    warning("iterative_forward (data): max_iter (", max_iter, ") reached ",
            "before convergence; inspect $history.", call. = FALSE)
  }
  anchor <- setdiff(linkable, flagged)
  list(
    rmsd = history[[length(history)]]$rmsd,
    rmsd_parameter = NULL, rmsd_data = history[[length(history)]]$rmsd,
    flagged = flagged, anchor = anchor,
    flagged_by_cutoff = stats::setNames(list(flagged), as.character(cut1)),
    linkable = linkable, history = history, converged = converged
  )
}
