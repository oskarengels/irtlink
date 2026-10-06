# R/dependent_vcov.R
# joint item parameter covariance for dependent person samples,
# method chooses the stacked sandwich or the person bootstrap

dependent_vcov <- function(x, ids, method = c("stacked", "bootstrap"),
                           ...) {
  method <- match.arg(method)
  if (inherits(x, "irtlink_calib") && x$model %in% c("GPCM", "PCM")) {
    stop("dependent_vcov is not yet available for polytomous models.",
         call. = FALSE)
  }
  if (inherits(x, "irtlink_calib") && !is.null(x$pweights)) {
    stop("dependent_vcov is not yet available for calibrations with ",
         "person weights.", call. = FALSE)
  }
  if (method == "stacked") {
    return(stacked_vcov(x, ids, ...))
  }
  bootstrap_vcov(x, ids, ...)
}
