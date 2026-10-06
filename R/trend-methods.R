# R/trend-methods.R

print.irtlink_trend <- function(x, ...) {
  link_lab <- if (!is.na(x$link)) paste0(", link = ", x$link) else ""
  cat("irtlink trend (", x$calibration, " / ", x$approach, link_lab, ")\n",
      sep = "")
  cat("Groups:", nrow(x$trend), "\n")
  if (!is.null(x$dif))
    cat("DIF purification: ", length(x$dif$flagged), " item(s) flagged\n",
        sep = "")
  cols <- c("group", "mu", "sigma")
  if ("le" %in% names(x$trend) && any(!is.na(x$trend$le))) {
    cat("Linking error: available (", sum(!is.na(x$trend$le)),
        " group estimate(s))\n", sep = "")
    cols <- c(cols, "le")
  }
  print(x$trend[cols], row.names = FALSE)
  invisible(x)
}

summary.irtlink_trend <- function(object, ...) {
  print(object)
  if (!is.null(object$dif) && length(object$dif$flagged))
    cat("\nFlagged items:", paste(object$dif$flagged, collapse = ", "), "\n")
  invisible(object)
}

coef.irtlink_trend <- function(object, ...) object$trend

plot.irtlink_trend <- function(x, ...) {
  tr <- x$trend
  do.call(graphics::plot, utils::modifyList(
    list(x = tr$group, y = tr$mu, type = "b", pch = 19,
         xlab = "group", ylab = expression(mu),
         ylim = range(c(tr$mu, 0), na.rm = TRUE),
         main = paste0("Trend (", x$calibration, "/", x$approach, ")")),
    list(...)))
  graphics::abline(h = 0, col = "grey70", lty = 3)
  invisible(x)
}
