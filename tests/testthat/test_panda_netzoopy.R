# runPANDA() gives the networks of netZooPy's Panda (CPU, double precision).
# References: fixtures/netzoopy, made by make_reference.py with netZooPy.

nzFile <- function(case, file) test_path("fixtures", "netzoopy", case, file)
nzMat <- function(case, file) {
  df <- utils::read.csv(nzFile(case, file), check.names = FALSE, stringsAsFactors = FALSE)
  m <- as.matrix(df[, -1, drop = FALSE])
  rownames(m) <- df[[1]]
  m
}
nzMeta <- function(case) {
  txt <- paste(readLines(nzFile(case, "meta.json"), warn = FALSE), collapse = "")
  list(mode = sub('.*"mode": *"([^"]*)".*', "\\1", txt),
       alpha = as.numeric(sub('.*"alpha": *([0-9.]+).*', "\\1", txt)))
}
# Same values for the same TF/gene names (row and column order may differ:
# MANTICORE sorts names with sort(), Python by code point).
expectNetZooPy <- function(x, ref, tol) {
  expect_setequal(rownames(x), rownames(ref))
  expect_setequal(colnames(x), colnames(ref))
  x <- x[rownames(ref), colnames(ref)]
  expect_lt(max(abs(x - ref)) / max(abs(ref)), tol)
}

test_that("runPANDA equals netZooPy's Panda", {
  for (case in c("basic", "union", "union_intersect", "noppi", "noexpr", "large")) {
    meta <- nzMeta(case)
    expr <- if (file.exists(nzFile(case, "expr.csv"))) nzMat(case, "expr.csv")
    motif <- utils::read.csv(nzFile(case, "motif.csv"), stringsAsFactors = FALSE)
    ppi <- if (file.exists(nzFile(case, "ppi.csv")))
      utils::read.csv(nzFile(case, "ppi.csv"), stringsAsFactors = FALSE)
    full <- file.exists(nzFile(case, "coexpression.csv"))
    r <- suppressMessages(MANTICORE:::runPANDA(
      motif = motif, expr = expr, ppi = ppi, alpha = meta$alpha, mode = meta$mode,
      output = c("regNet", "coregNet", "coopNet"), progress = FALSE))
    expectNetZooPy(r$regNet, nzMat(case, "regulatory.csv"), 1e-12)
    if (full) {
      expectNetZooPy(r$coopNet, nzMat(case, "cooperative.csv"), 1e-12)
      expectNetZooPy(r$coregNet, nzMat(case, "coexpression.csv"), 1e-12)
    }
  }
})

test_that("PPI and motif are built as netZooPy builds them", {
  motif <- data.frame(tf = c("b", "a", "a", "a"), gene = c("y", "x", "z", "x"), w = c(1, 2, 3, 4))
  ppi <- data.frame(a = c("a", "b"), b = c("b", "a"), w = c(0.9, 0.3))
  expr <- matrix(1:12, 3, 4, dimnames = list(c("z", "x", "q"), NULL))
  i <- MANTICORE:::.pandaInputs(motif, expr, ppi, "intersection")
  expect_identical(i$regulatoryNetwork["a", "x"], 4)   # last weight wins
  expect_identical(i$tfCoopNetwork["a", "b"], 0.3)     # second pass wins
  expect_identical(i$tfCoopNetwork["b", "a"], 0.9)
  expect_identical(unname(diag(i$tfCoopNetwork)), c(1, 1))
  u <- MANTICORE:::.pandaInputs(motif, expr, ppi, "union")
  expect_identical(u$gene.names, c("q", "x", "y", "z"))
  expect_identical(unname(u$expr["z", ]), c(1, 4, 7, 10))   # rows matched by name
  expect_identical(unname(u$expr["y", ]), rep(0, 4))
  self <- MANTICORE:::.pandaInputs(motif, expr, rbind(ppi, data.frame(a = "a", b = "a", w = 0.5)),
                                  "intersection")
  expect_identical(self$tfCoopNetwork["a", "a"], 0.5)
})

test_that("the diagonal update uses the population SD, as netZooPy", {
  m <- matrix(c(9, 1, 2, 3, 9, 5, 6, 7, 9), 3)
  d <- MANTICORE:::update.diagonal(m, 3, 0.1, 2)
  sdPop <- function(x) sqrt(mean((x - mean(x))^2))
  expect_equal(diag(d), c(sdPop(c(1, 2)), sdPop(c(3, 5)), sdPop(c(6, 7))) * 3 * exp(0.4))
})
