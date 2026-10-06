print.irtlink_link <- function(x, ...) {
  cat("irtlink ", x$approach, " linking (", x$method, ")\n", sep = "")
  if (!is.null(x$steps)) {
    cat("Steps (n common items):",
        paste(sprintf("%d->%d (%d)", x$steps$from, x$steps$to,
                      x$steps$n_common), collapse = ", "), "\n")
    not_conv <- which(is.na(x$steps$converged) | !x$steps$converged)
    if (length(not_conv) > 0) {
      cat("WARNING: step(s) did not converge:",
          paste(not_conv, collapse = ", "), "\n")
    }
  } else {
    cat("Groups:", nrow(x$trend),
        "| Items:", length(unique(x$ipars$item)), "\n")
    if (any(is.na(x$converged) | !x$converged)) {
      cat("WARNING: joint estimation did not (verifiably) converge\n")
    }
  }
  tr <- x$trend
  tr$mu <- round(tr$mu, 3)
  tr$sigma <- round(tr$sigma, 3)
  print(tr, row.names = FALSE)
  invisible(x)
}

summary.irtlink_link <- function(object, ...) {
  print(object)
  if (!is.null(object$steps)) {
    cat("\nStep details:\n")
    print(object$steps, row.names = FALSE)
  }
  invisible(object)
}

coef.irtlink_link <- function(object, ...) {
  object$trend
}

plot.irtlink_link <- function(x, pars = c("b", "a"), ...) {
  pars <- match.arg(pars)
  if (is.null(x$ipars)) {
    stop("This irtlink_link object carries no item parameters; ",
         "re-run the linking to enable plotting.", call. = FALSE)
  }
  tr <- x$trend
  n_steps <- nrow(tr) - 1
  old <- graphics::par(mfrow = grDevices::n2mfrow(n_steps))
  on.exit(graphics::par(old), add = TRUE, after = FALSE)
  for (t in seq_len(n_steps)) {
    # joint links may connect neighboring groups only indirectly
    pair <- tryCatch(ipars_pair(x$ipars, tr$group[t], tr$group[t + 1]),
                     error = function(e) NULL)
    if (is.null(pair)) {
      graphics::plot.new()
      graphics::title(main = sprintf("%d -> %d: no common items",
                                     tr$group[t], tr$group[t + 1]))
      next
    }
    # trend-implied step transformation for chain and joint results
    sigma_step <- tr$sigma[t + 1] / tr$sigma[t]
    mu_step <- (tr$mu[t + 1] - tr$mu[t]) / tr$sigma[t]
    if (pars == "b") {
      xv <- pair$b_from
      yv <- mu_step + sigma_step * pair$b_to
      lab <- "difficulty b"
    } else {
      xv <- pair$a_from
      yv <- pair$a_to / sigma_step
      lab <- "discrimination a"
    }
    do.call(graphics::plot, c(
      list(x = xv, y = yv),
      utils::modifyList(
        list(xlab = sprintf("group %d", tr$group[t]),
             ylab = sprintf("group %d (transformed)", tr$group[t + 1]),
             main = sprintf("%d -> %d: %s", tr$group[t],
                            tr$group[t + 1], lab)),
        list(...)
      )
    ))
    graphics::abline(0, 1, lty = 2)
  }
  invisible(x)
}
