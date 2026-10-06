# R/em-parmap.R
# parameter map for the internal EM engine, one row per
# (item, group, parameter) cell

em_parmap_separate <- function(items, model) {
  reg <- em_model_registry(model)
  J <- length(items)
  np <- length(reg$parnames)
  pm <- data.frame(
    item = rep(items, each = np),
    group = 1L,
    parname = rep(reg$parnames, J),
    parindex = NA_integer_,
    est = rep(unname(reg$est), J),
    value = rep(unname(reg$start), J),
    lower = rep(unname(reg$lower), J),
    upper = rep(unname(reg$upper), J),
    penalty_group = NA_character_,
    stringsAsFactors = FALSE
  )
  pm$parindex[pm$est] <- seq_len(sum(pm$est))
  pm
}

# concurrent map, a and nu tied across groups, DIF offsets free
# only at later administrations of freed items
em_parmap_concurrent <- function(md, model, free_items = NULL,
                                 free_slope_items = NULL) {
  reg <- em_model_registry(model)
  freed_any <- union(free_items, free_slope_items)
  if (length(freed_any) > 0) {
    unknown <- setdiff(freed_any, md$pool)
    if (length(unknown) > 0)
      warning("free_items not found in data: ",
              paste(unknown, collapse = ", "), call. = FALSE)
  }
  rows <- list()
  next_index <- 0L
  for (jj in seq_along(md$pool)) {
    it <- md$pool[jj]
    tps <- which(vapply(md$groups, function(w) jj %in% w$idx, logical(1)))
    idx_a <- if (reg$est[["a"]]) next_index + 1L else NA_integer_
    idx_nu <- next_index + 1L + reg$est[["a"]]
    next_index <- next_index + reg$est[["a"]] + 1L
    for (t in tps) {
      # the first administration anchors (a, nu)
      g_free <- (it %in% free_items) && t != tps[1L]
      if (g_free) next_index <- next_index + 1L
      idx_g <- if (g_free) next_index else NA_integer_
      h_free <- (it %in% free_slope_items) && t != tps[1L]
      if (h_free) next_index <- next_index + 1L
      idx_h <- if (h_free) next_index else NA_integer_
      rows[[length(rows) + 1L]] <- data.frame(
        item = it, group = t,
        parname = c("a", "nu", "g", "h"),
        parindex = c(idx_a, idx_nu, idx_g, idx_h),
        est = c(reg$est[["a"]], TRUE, g_free, h_free),
        value = c(unname(reg$start[["a"]]), unname(reg$start[["nu"]]),
                  0, 0),
        lower = c(unname(reg$lower), -10, -10),
        upper = c(unname(reg$upper), 10, 10),
        penalty_group = NA_character_,
        stringsAsFactors = FALSE
      )
    }
  }
  do.call(rbind, rows)
}
