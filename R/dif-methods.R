# R/dif-methods.R

print.irtlink_dif <- function(x, ...) {
  if (identical(x$method, "lrt")) {
    adj <- if (!is.null(x$alpha_adjust) && x$alpha_adjust == "bonferroni")
      ", Bonferroni" else ""
    cat("irtlink DIF detection (LRT, dif_type = ", x$dif_type, adj, ")\n",
        sep = "")
    cat("alpha:", x$alpha, "| min_anchor:", x$min_anchor, "\n")
    fl_str <- if (length(x$flagged)) {
      paste0(" (", paste(x$flagged, collapse = ", "), ")")
    } else {
      ""
    }
    cat("Flagged: ", length(x$flagged), fl_str, "\n", sep = "")
    cat("Anchor:", length(x$anchor), "of", length(x$linkable),
        "linkable items\n")
    if (!is.null(x$lrt)) {
      pcol <- if (!is.null(x$alpha_adjust) && x$alpha_adjust == "bonferroni" &&
                  "p_adj" %in% names(x$lrt)) "p_adj" else "p"
      ord <- order(x$lrt[[pcol]])[seq_len(min(5L, nrow(x$lrt)))]
      cat("Smallest ", if (pcol == "p_adj") "adjusted " else "", "p-values:\n",
          sep = "")
      print(x$lrt[ord, c("item", "chisq", "df", pcol)], row.names = FALSE)
    }
    return(invisible(x))
  }
  cat("irtlink DIF detection (", toupper(x$method), ", ",
      x$rmsd_method, ", ", x$procedure, ")\n", sep = "")
  cut_lab <- if (x$cutoff_type == "fixed") {
    paste0("fixed cutoff ", paste(x$cutoff, collapse = "/"))
  } else {
    paste0("data-driven tau ", x$tau)
  }
  cat("Criterion:", cut_lab, "| min_anchor:", x$min_anchor, "\n")
  fl_str <- if (length(x$flagged)) {
    paste0(" (", paste(x$flagged, collapse = ", "), ")")
  } else {
    ""
  }
  cat("Flagged: ", length(x$flagged), fl_str, "\n", sep = "")
  cat("Anchor:", length(x$anchor), "of", length(x$linkable),
      "linkable items\n")
  invisible(x)
}

summary.irtlink_dif <- function(object, ...) {
  print(object)
  if (!is.null(object$flagged_by_cutoff) &&
      length(object$flagged_by_cutoff) > 1) {
    cat("\nFlagging by cutoff:\n")
    for (nm in names(object$flagged_by_cutoff)) {
      fl <- object$flagged_by_cutoff[[nm]]
      cat("  ", nm, ": ",
          if (length(fl)) paste(fl, collapse = ", ") else "(none)", "\n",
          sep = "")
    }
  }
  if (!is.null(object$rmsd)) {
    cat("\nTop RMSD items:\n")
    tab <- object$rmsd[order(object$rmsd$max, decreasing = TRUE), ]
    print(utils::head(tab[c("item", "max")], 10), row.names = FALSE)
  } else if (!is.null(object$lrt)) {
    cat("\nTop LRT items (by chi-square):\n")
    tab <- object$lrt[order(object$lrt$chisq, decreasing = TRUE), ]
    print(utils::head(tab[c("item", "chisq", "df", "p")], 10),
          row.names = FALSE)
  }
  invisible(object)
}

coef.irtlink_dif <- function(object, ...) {
  if (!is.null(object$rmsd)) object$rmsd else object$lrt
}

plot.irtlink_dif <- function(x, ...) {
  if (identical(x$method, "lrt")) {
    tab <- x$lrt[order(x$lrt$chisq, decreasing = TRUE), ]
    if (nrow(tab) == 0L) {
      return(invisible(x))
    }
    cols <- ifelse(tab$item %in% x$flagged, "red", "black")
    do.call(graphics::plot, utils::modifyList(
      list(x = seq_len(nrow(tab)), y = tab$chisq, col = cols, pch = 19,
           xlab = "item (sorted by LRT chi-square)",
           ylab = "LRT chi-square", main = "DIF detection (LRT)"),
      list(...)
    ))
    # the df can differ across items administered in different
    # numbers of groups
    if (length(unique(tab$df)) == 1L) {
      graphics::abline(h = stats::qchisq(1 - x$alpha, df = tab$df[1]),
                       lty = 2)
    } else {
      crit <- stats::qchisq(1 - x$alpha, df = tab$df)
      graphics::segments(seq_len(nrow(tab)) - 0.4, crit,
                         seq_len(nrow(tab)) + 0.4, crit, lty = 2)
    }
    return(invisible(x))
  }
  tab <- x$rmsd[order(x$rmsd$max, decreasing = TRUE), ]
  if (nrow(tab) == 0L) {
    return(invisible(x))
  }
  flagged <- tab$item %in% x$flagged
  cols <- ifelse(flagged, "red", "black")
  do.call(graphics::plot, utils::modifyList(
    list(x = seq_len(nrow(tab)), y = tab$max, col = cols, pch = 19,
         xlab = "item (sorted by RMSD)", ylab = "max RMSD across groups",
         main = "DIF detection"),
    list(...)
  ))
  if (x$cutoff_type == "fixed") {
    graphics::abline(h = x$cutoff[1], lty = 2)
  }
  invisible(x)
}
