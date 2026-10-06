# R/dif-flag.R

# keep at least min_anchor anchor items when flagging DIF
dif_tiebreak <- function(stat, above, min_anchor, total) {
  n_remaining <- total - length(above)
  if (n_remaining >= min_anchor) {
    return(above)
  }
  max_allowed <- total - min_anchor
  if (max_allowed <= 0) {
    return(character(0))
  }
  # ties broken by row order
  ordered <- names(sort(stat[above], decreasing = TRUE))
  ordered[seq_len(min(max_allowed, length(ordered)))]
}

# fixed-cutoff flagging on the per-item RMSD maximum
dif_flag_fixed <- function(tab, cutoff, min_anchor) {
  stat <- stats::setNames(tab$max, tab$item)
  above <- tab$item[tab$max > cutoff]
  dif_tiebreak(stat, above, min_anchor, total = nrow(tab))
}

# robust z score against the own median, zero MAD gives z = 0
robust_z <- function(x) {
  med <- stats::median(x, na.rm = TRUE)
  madv <- stats::median(abs(x - med), na.rm = TRUE) * 1.4826
  if (!is.finite(madv) || madv == 0) {
    return(ifelse(is.na(x), NA_real_, 0))
  }
  abs(x - med) / madv
}

# data-driven flagging with the robust z score
# (von Davier & Bezirhan, 2023)
dif_flag_mad <- function(tab, tau, min_anchor,
                         group_cols = grep("^g[0-9]+$", names(tab),
                                          value = TRUE)) {
  z <- vapply(group_cols, function(cc) robust_z(tab[[cc]]),
              numeric(nrow(tab)))
  z <- matrix(z, nrow = nrow(tab))
  z_max <- apply(z, 1, max, na.rm = TRUE)
  stat <- stats::setNames(z_max, tab$item)
  above <- tab$item[z_max > tau]
  dif_tiebreak(stat, above, min_anchor, total = nrow(tab))
}
