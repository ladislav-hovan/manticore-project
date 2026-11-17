# As netZooPy's update_diagonal: the diagonal becomes the population SD
# (numpy.nanstd, ddof = 0) of each column's other values, times the number of
# nodes, times exp(2 * alpha * step).
update.diagonal <- function(
  diagMat,
  num,
  alpha,
  step
) {
  diag(diagMat) <- NaN
  diagstd <- .colStats(diagMat, na.rm = TRUE)$sd
  diag(diagMat) <- diagstd * num * exp(2 * alpha * step)

  return(diagMat)
}