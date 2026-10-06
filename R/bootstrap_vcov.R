# R/bootstrap_vcov.R
# Nonparametric person bootstrap of the joint item parameter covariance
# for longitudinal designs in which the person samples of the groups
# overlap. Estimates the same matrix as stacked_vcov().

bootstrap_vcov <- function(x, ids, B = 200, seed = NULL, model = "2PL",
                           cores = 1L, draws = NULL, ...) {
  if (!is.list(ids))
    stop("`ids` must be a list with one identifier vector per group.",
         call. = FALSE)
  if (is.null(draws)) {
    if (!is.list(x) || length(x) < 2L)
      stop("`x` must be a list with one response matrix per group.",
           call. = FALSE)
    Tn <- length(x)
    if (!is.list(ids) || length(ids) != Tn)
      stop("`ids` must be a list with one identifier vector per group.",
           call. = FALSE)
    ids <- lapply(ids, as.character)
    resp_list <- lapply(x, as.matrix)
    for (t in seq_len(Tn)) {
      if (length(ids[[t]]) != nrow(resp_list[[t]]))
        stop("`ids[[", t, "]]` must have one identifier per row of the ",
             "group-", t, " response data (", nrow(resp_list[[t]]),
             " rows).", call. = FALSE)
      if (anyDuplicated(ids[[t]]))
        stop("`ids[[", t, "]]` contains duplicated identifiers.",
             call. = FALSE)
    }
    B <- as.integer(B)
    if (is.na(B) || B < 2L)
      stop("`B` must be at least 2.", call. = FALSE)
    cores <- max(1L, as.integer(cores))
    if (cores > 1L && is.null(seed))
      stop("Reproducible parallel draws need `seed`; set `seed` or use ",
           "cores = 1.", call. = FALSE)

    union_ids <- unique(unlist(ids))
    N_union <- length(union_ids)
    row_of <- lapply(ids, function(v) stats::setNames(seq_along(v), v))

    boot_one <- function(b) {
      if (!is.null(seed)) set.seed(seed + b)
      draw <- sample(union_ids, N_union, replace = TRUE)
      resp_b <- lapply(seq_along(resp_list), function(k) {
        rows <- row_of[[k]][draw]
        rows <- rows[!is.na(rows)]
        X <- resp_list[[k]][rows, , drop = FALSE]
        rownames(X) <- NULL
        v <- apply(X, 2, function(col) stats::var(col, na.rm = TRUE))
        drop <- is.na(v) | v == 0
        list(X = X[, !drop, drop = FALSE], dropped = colnames(X)[drop])
      })
      dropped <- unlist(lapply(resp_b, `[[`, "dropped"))
      cal_b <- suppressWarnings(calibrate(lapply(resp_b, `[[`, "X"),
                                          model = model,
                                          calibration = "separate", ...))
      ip <- cal_b$ipars
      data.frame(rep = b, group = ip$group, item = ip$item,
                 a = ip$a, b = ip$b,
                 converged = paste(cal_b$converged, collapse = ","),
                 dropped = paste(dropped, collapse = ";"),
                 stringsAsFactors = FALSE)
    }

    if (cores == 1L) {
      draw_list <- lapply(seq_len(B), boot_one)
    } else {
      cl <- parallel::makeCluster(cores)
      on.exit(parallel::stopCluster(cl), add = TRUE)
      parallel::clusterCall(cl, function(paths) {
        .libPaths(paths)
        loadNamespace("irtlink")
        NULL
      }, .libPaths())
      draw_list <- parallel::parLapplyLB(cl, seq_len(B), boot_one)
    }
    draws <- do.call(rbind, draw_list)
  } else {
    need <- c("rep", "group", "item", "a", "b", "converged", "dropped")
    if (!is.data.frame(draws) || !all(need %in% names(draws)))
      stop("`draws` must be a data frame with the columns ",
           paste(need, collapse = ", "), ".", call. = FALSE)
    B <- length(unique(draws$rep))
  }
  overlap <- vcov_overlap_matrix(lapply(ids, as.character))

  nonconv <- unique(draws$rep[
    vapply(strsplit(draws$converged, ","),
           function(v) any(v != "TRUE"), logical(1))])
  with_drop <- unique(draws$rep[!is.na(draws$dropped) &
                                  draws$dropped != ""])
  excluded <- sort(union(nonconv, with_drop))
  used <- draws[!draws$rep %in% excluded, , drop = FALSE]
  if (length(unique(used$rep)) < 2L)
    stop("Fewer than two usable replicates; increase `B`.", call. = FALSE)

  groups <- sort(unique(used$group))
  pars <- unlist(lapply(groups, function(t) {
    its <- unique(used$item[used$group == t])
    as.vector(t(cbind(paste0("g", t, ":", its, ":a"),
                      paste0("g", t, ":", its, ":b"))))
  }))
  reps <- sort(unique(used$rep))
  G <- matrix(NA_real_, length(reps), length(pars),
              dimnames = list(reps, pars))
  pa <- paste0("g", used$group, ":", used$item, ":a")
  pb <- paste0("g", used$group, ":", used$item, ":b")
  G[cbind(match(used$rep, reps), match(pa, pars))] <- used$a
  G[cbind(match(used$rep, reps), match(pb, pars))] <- used$b
  complete <- stats::complete.cases(G)
  G <- G[complete, , drop = FALSE]
  V <- stats::cov(G)
  V <- (V + t(V)) / 2

  items_all <- sort(unique(used$item))
  vcov_items <- lapply(items_all, function(it) {
    nm <- pars[grepl(paste0(":", it, ":"), pars, fixed = TRUE)]
    V[nm, nm, drop = FALSE]
  })
  names(vcov_items) <- items_all

  structure(
    list(vcov = V, vcov_items = vcov_items, draws = draws,
         B = B, B_used = nrow(G), excluded = excluded,
         overlap = overlap, n_groups = length(groups),
         call = match.call()),
    class = "irtlink_bootstrap_vcov"
  )
}

print.irtlink_bootstrap_vcov <- function(x, ...) {
  cat("Person-bootstrap covariance of the item parameters (",
      x$n_groups, " groups, ", nrow(x$vcov), " parameters)\n", sep = "")
  cat("Replicates used: ", x$B_used, " of ", x$B, sep = "")
  if (length(x$excluded)) {
    cat(" (excluded: ", paste(x$excluded, collapse = ", "), ")", sep = "")
  }
  cat("\n")
  if (!is.null(x$overlap)) {
    cat("Shared persons per group pair:\n")
    print(x$overlap)
  }
  invisible(x)
}
