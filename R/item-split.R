# R/item-split.R
# split flagged items into one column per group for partial
# invariance fits

# split `split_items` in a stacked response matrix X
split_items_resp <- function(X, group, split_items, n_groups) {
  keep <- setdiff(colnames(X), split_items)
  out <- X[, keep, drop = FALSE]
  if (length(split_items) == 0L) return(out)
  add <- lapply(split_items, function(it) {
    cols <- vapply(seq_len(n_groups), function(g) {
      v <- rep(NA_real_, nrow(X))
      sel <- group == g
      v[sel] <- X[sel, it]
      v
    }, numeric(nrow(X)))
    colnames(cols) <- paste0(it, "__G", seq_len(n_groups))
    cols[, vapply(seq_len(ncol(cols)), function(j) any(!is.na(cols[, j])),
                  logical(1)), drop = FALSE]
  })
  cbind(out, do.call(cbind, add))
}

# raw responses from one fitted model
calib_raw_responses <- function(m) {
  r <- if (!is.null(m$resp) && (is.matrix(m$resp) || is.data.frame(m$resp))) {
    m$resp
  } else if (!is.null(m$dat) && (is.matrix(m$dat) || is.data.frame(m$dat))) {
    m$dat
  } else {
    stop("Could not locate the raw responses in the fitted model ",
         "(neither $resp (matrix) nor $dat (matrix) found).", call. = FALSE)
  }
  as.data.frame(r)
}

# per-group response data from the stored models
responses_per_group <- function(calib) {
  if (!identical(calib$calibration, "separate"))
    stop("This detection method requires separate per-group calibration ",
         "(calibration = \"separate\").", call. = FALSE)
  if (is.null(calib$models))
    stop("Needs the fitted models; re-run calibrate() with keep_models = TRUE.",
         call. = FALSE)
  lapply(calib$models, calib_raw_responses)
}
