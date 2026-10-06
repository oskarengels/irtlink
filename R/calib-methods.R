# R/calib-methods.R

print.irtlink_calib <- function(x, ...) {
  cat("irtlink calibration (", x$model, ", ", x$calibration, ")\n", sep = "")
  np <- if (!is.null(x$N)) paste(x$N, collapse = ", ") else "NA"
  cat("Groups:", x$n_groups, "| Persons per group:", np, "\n")
  # ipars can be NULL for composed concurrent fits
  if (!is.null(x$ipars)) {
    items_per_group <- table(x$ipars$group)
    cat("Items per group:",
        paste(names(items_per_group), as.integer(items_per_group),
              sep = ":", collapse = "  "), "\n")
  }
  if (!is.null(x$trend)) {
    cat("Trend (mu, sigma) per group:\n")
    cat(paste(sprintf("  g%d: mu=%.3f sigma=%.3f",
                      x$trend$group, x$trend$mu, x$trend$sigma),
              collapse = "\n"), "\n")
    if (!is.null(x$free_items) && length(x$free_items) > 0)
      cat("Freed (partial-invariance) items:",
          paste(x$free_items, collapse = ", "), "\n")
    if (!is.null(x$eps))
      cat("SBIC eps (selected):", x$eps, "\n")
  }
  # NA convergence (unknown) is reported alongside non-convergence
  not_conv <- which(is.na(x$converged) | !x$converged)
  if (length(not_conv) > 0) {
    labels <- if (!is.null(names(x$converged))) names(x$converged)[not_conv] else not_conv
    cat("WARNING: calibration did not converge for:",
        paste(labels, collapse = ", "), "\n")
  }
  invisible(x)
}

summary.irtlink_calib <- function(object, ...) {
  print(object)
  if (is.null(object$ipars)) return(invisible(object))
  if (anyNA(object$ipars[c("a", "b")])) {
    warning("ipars contains NA - parameter ranges exclude affected items.",
            call. = FALSE)
  }
  cat("\nItem parameter ranges per group:\n")
  rng <- aggregate(cbind(a, b) ~ group, data = object$ipars,
                   FUN = function(v) c(min = min(v), max = max(v)))
  print(rng)
  invisible(object)
}

coef.irtlink_calib <- function(object,
                               type = c("difficulty", "intercept", "full"),
                               ...) {
  type <- match.arg(type)
  if (type == "difficulty") {
    return(object$ipars)
  }
  if (is.null(object$item_parameters)) {
    stop("This calibration object does not contain `item_parameters`.",
         call. = FALSE)
  }
  if (type == "full") {
    return(object$item_parameters)
  }
  object$item_parameters[c("group", "item", "a", "d")]
}
