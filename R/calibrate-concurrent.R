# R/calibrate-concurrent.R
# multiple group (concurrent) calibration

# per-group item parameters of a multiple group fit, freed items
# carry their per-group offsets (a + h, nu + g)
concurrent_ipars <- function(fit, data_list) {
  out <- lapply(seq_along(data_list), function(t) {
    items_t <- colnames(as.data.frame(data_list[[t]]))
    idx <- match(items_t, fit$item)
    a_t <- fit$a[idx] + fit$h[idx, t]
    nu_t <- fit$nu[idx] + fit$g[idx, t]
    data.frame(group = t, item = items_t, a = a_t, b = nu_t / a_t,
               stringsAsFactors = FALSE)
  })
  do.call(rbind, out)
}
