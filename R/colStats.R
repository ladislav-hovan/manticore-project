# Means and population SDs (two-pass, as numpy and scipy compute them) of the
# columns or rows of a matrix, without full-size temporary matrices.
.colStats <- function(X, na.rm = FALSE) {
  mu <- colMeans(X, na.rm = na.rm)
  st <- vapply(seq_len(ncol(X)), function(j) {
    x <- X[, j] - mu[j]
    c(sum(x * x, na.rm = na.rm), if (na.rm) sum(!is.na(x)) else length(x))
  }, numeric(2))
  list(mean = mu, sd = sqrt(st[1, ] / st[2, ]))
}

.rowStats <- function(X, block = 256L) {
  mu <- rowMeans(X)
  ss <- numeric(nrow(X))
  for (j0 in seq(1L, ncol(X), by = block)) {
    J <- j0:min(ncol(X), j0 + block - 1L)
    ss <- ss + rowSums((X[, J, drop = FALSE] - mu)^2)
  }
  list(mean = mu, sd = sqrt(ss / ncol(X)))
}
