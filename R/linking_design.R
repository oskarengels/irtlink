linking_design <- function(data_list) {
  if (!is.list(data_list) || length(data_list) < 2) {
    stop("`data_list` must be a list of at least two data sets ",
         "(one per group).", call. = FALSE)
  }
  nm <- lapply(seq_along(data_list), function(t) {
    cn <- colnames(data_list[[t]])
    if (is.null(cn)) {
      stop("Data set for group ", t, " has no column names; ",
           "items must be named.", call. = FALSE)
    }
    cn
  })
  pool <- sort(unique(unlist(nm)))
  mat <- vapply(nm, function(items) pool %in% items,
                logical(length(pool)))
  if (!is.matrix(mat)) {
    # vapply collapses to a vector when the pool has a single item
    dim(mat) <- c(length(pool), length(nm))
  }
  rownames(mat) <- pool
  colnames(mat) <- paste0("group", seq_along(data_list))
  successive <- vapply(seq_len(ncol(mat) - 1), function(t) {
    sum(mat[, t] & mat[, t + 1])
  }, integer(1))
  structure(
    list(matrix = mat, successive_common = successive,
         n_groups = length(data_list)),
    class = "irtlink_design"
  )
}

print.irtlink_design <- function(x, ...) {
  cat("irtlink linking design:", x$n_groups, "groups,",
      nrow(x$matrix), "items in pool\n")
  cat("Common items between successive groups:",
      paste(x$successive_common, collapse = ", "), "\n")
  if (nrow(x$matrix) <= 60) {
    m <- ifelse(x$matrix, "x", ".")
    print(as.data.frame(m))
  } else {
    cat("(matrix suppressed for pools > 60 items)\n")
  }
  invisible(x)
}
