# R/pinv.R
# Moore-Penrose pseudoinverse via the singular value decomposition
pinv <- function(A, tol = sqrt(.Machine$double.eps)) {
  s <- svd(A)
  keep <- s$d > tol * s$d[1]
  if (!any(keep)) return(matrix(0, ncol(A), nrow(A)))
  V <- sweep(s$v[, keep, drop = FALSE], 2L, s$d[keep], "/")
  tcrossprod(V, s$u[, keep, drop = FALSE])
}
