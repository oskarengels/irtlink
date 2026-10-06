# internal long-format item parameter representation, one row per
# item and group, all linking functions operate on this format

as_ipars <- function(x) {
  if (inherits(x, "irtlink_calib")) {
    ip <- x$ipars
    if (!is.null(ip) && is.null(ip$c)) ip$c <- 0
    return(ip)
  }
  if (!is.data.frame(x)) {
    stop("`x` must be an `irtlink_calib` object or a data.frame ",
         "with columns group, item, a, b.", call. = FALSE)
  }
  required <- c("group", "item", "a", "b")
  missing_cols <- setdiff(required, names(x))
  if (length(missing_cols) > 0) {
    stop("Missing columns in item parameter data: ",
         paste(missing_cols, collapse = ", "), call. = FALSE)
  }
  has_c <- "c" %in% names(x)
  tau_cols <- grep("^tau[0-9]+$", names(x), value = TRUE)
  x <- x[c(required, if (has_c) "c", tau_cols)]
  if (!has_c) x$c <- 0
  if (!is.numeric(x$c) || anyNA(x$c) || any(x$c < 0) || any(x$c >= 1)) {
    stop("Column `c` (guessing) must be numeric in [0, 1).", call. = FALSE)
  }
  for (tc in tau_cols) {
    if (!is.numeric(x[[tc]])) {
      stop("Threshold column `", tc, "` must be numeric.", call. = FALSE)
    }
  }
  x$group <- as.integer(x$group)
  x$item <- as.character(x$item)
  if (!is.numeric(x$a) || !is.numeric(x$b)) {
    stop("Columns `a` and `b` must be numeric.", call. = FALSE)
  }
  if (anyNA(x$group)) {
    stop("Column `group` must not contain NA.", call. = FALSE)
  }
  if (anyNA(x$a) || anyNA(x$b)) {
    stop("Columns `a` and `b` must not contain NA.", call. = FALSE)
  }
  groups <- sort(unique(x$group))
  if (!identical(groups, seq_along(groups))) {
    stop("`group` must contain consecutive integers starting at 1.",
         call. = FALSE)
  }
  x <- x[order(x$group, x$item), ]
  if (anyDuplicated(x[c("group", "item")])) {
    stop("Duplicate (group, item) combinations found in item parameter data.",
         call. = FALSE)
  }
  rownames(x) <- NULL
  x
}

# TRUE when the item parameters carry polytomous thresholds
ipars_poly <- function(ipars) {
  tc <- grep("^tau[0-9]+$", names(ipars))
  length(tc) > 0 && any(!is.na(as.matrix(ipars[tc])))
}

# wide format for one linking step, common items only
ipars_pair <- function(ipars, group_from, group_to) {
  from <- ipars[ipars$group == group_from, ]
  to <- ipars[ipars$group == group_to, ]
  common <- intersect(from$item, to$item)
  if (length(common) == 0) {
    stop(sprintf("No common items between group %d and group %d.",
                 group_from, group_to), call. = FALSE)
  }
  from <- from[match(common, from$item), ]
  to <- to[match(common, to$item), ]
  pair <- data.frame(item = common,
                     a_from = from$a, b_from = from$b,
                     a_to = to$a, b_to = to$b,
                     c_from = if (is.null(from$c)) 0 else from$c,
                     c_to = if (is.null(to$c)) 0 else to$c)
  tau_cols <- grep("^tau[0-9]+$", names(ipars), value = TRUE)
  for (tc in tau_cols) {
    pair[[paste0(tc, "_from")]] <- from[[tc]]
    pair[[paste0(tc, "_to")]] <- to[[tc]]
  }
  pair
}
