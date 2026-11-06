# As netZooPy's normalize_network: (column z-score + row z-score) / sqrt(2).
# A z-score is NaN where scipy's zscore() makes it NaN (constant rows or
# columns, SD <= eps * |mean|); there the z-score against the whole matrix
# (sample SD) is used instead.
normalizeNetwork <- function(X) {
  nr <- nrow(X)
  nc <- ncol(X)

  # Row and column stats (vectors only — no full-matrix broadcast)
  rs <- .rowStats(X)
  cs <- .colStats(X)
  mu1 <- rs$mean
  std1 <- rs$sd
  mu2 <- cs$mean
  std2 <- cs$sd
  eps <- .Machine$double.eps
  flat1 <- std1 <= abs(eps * mu1)
  flat2 <- std2 <= abs(eps * mu2)

  if (!any(flat1) && !any(flat2)) {
    # Fast path: accumulate Z1 + Z2 in a single matrix
    normMat <- (X - mu1) / std1
    for (j in seq_len(nc)) {
      normMat[, j] <- normMat[, j] + (X[, j] - mu2[j]) / std2[j]
    }
    normMat <- normMat / sqrt(2)
  } else {
    # Slow path: z-scores of constant rows/columns are NaN
    Z1 <- (X - mu1) / std1
    Z1[flat1, ] <- NaN

    Z2 <- X
    for (j in seq_len(nc)) {
      Z2[, j] <- if (flat2[j]) NaN else (X[, j] - mu2[j]) / std2[j]
    }

    normMat <- (Z1 + Z2) / sqrt(2)

    Z0 <- (X - mean(X)) / stats::sd(as.vector(X))
    f1 <- is.na(Z1)
    f2 <- is.na(Z2)

    normMat[f1] <- (Z2[f1] + Z0[f1]) / sqrt(2)
    normMat[f2] <- (Z1[f2] + Z0[f2]) / sqrt(2)
    normMat[f1 & f2] <- 2 * Z0[f1 & f2] / sqrt(2)
  }

  normMat
}
