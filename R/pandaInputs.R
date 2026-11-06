# Builds PANDA's starting matrices with the same rules as netZooPy's
# Panda.processData (modes "intersection" and "union"), so runPANDA() gives
# the networks netZooPy gives:
# * genes: in the motif and the expression ("intersection") or in either
#   ("union", with zero expression for genes that have none);
#   TFs: in the motif and the PPI, or in either. Without a PPI, the motif TFs;
#   without expression, the motif genes.
# * a (TF, gene) pair given twice in the motif keeps the last weight;
# * the PPI matrix starts as the identity and is filled in both directions:
#   first every (a, b) row sets [a, b], then every row sets [b, a]. A pair
#   given in both orders with different weights is therefore asymmetric, and
#   the diagonal is 1 unless a self-pair gives another weight.
# Names are sorted with sort(), as SCORPION always has.
.pandaInputs <- function(motif, expr, ppi, mode) {
  motifTfs <- unique(as.character(motif[, 1]))
  motifGenes <- unique(as.character(motif[, 2]))
  ppiTfs <- if (is.null(ppi)) motifTfs else unique(as.character(c(ppi[, 1], ppi[, 2])))
  exprGenes <- if (is.null(expr)) motifGenes else rownames(expr)
  if (mode == "union") {
    gene.names <- sort(union(motifGenes, exprGenes))
    tf.names <- sort(union(ppiTfs, motifTfs))
  } else {
    gene.names <- sort(intersect(exprGenes, motifGenes))
    tf.names <- sort(intersect(ppiTfs, motifTfs))
  }
  num.TFs <- length(tf.names)
  num.genes <- length(gene.names)

  if (!is.null(expr)) {
    if (anyDuplicated(rownames(expr))) {
      cli::cli_abort("The expression has repeated gene names. Average or remove them first.")
    }
    if (mode == "union") {
      expr1 <- matrix(0, num.genes, ncol(expr), dimnames = list(gene.names, colnames(expr)))
      hit <- match(rownames(expr), gene.names)
      expr1[hit, ] <- as.matrix(expr)
      expr <- expr1
    } else {
      expr <- expr[gene.names, , drop = FALSE]
    }
  }

  regulatoryNetwork <- matrix(0, num.TFs, num.genes, dimnames = list(tf.names, gene.names))
  Idx1 <- match(motif[, 1], tf.names)
  Idx2 <- match(motif[, 2], gene.names)
  ok <- !is.na(Idx1) & !is.na(Idx2)
  regulatoryNetwork[cbind(Idx1[ok], Idx2[ok])] <- as.numeric(motif[ok, 3])

  tfCoopNetwork <- diag(num.TFs)
  dimnames(tfCoopNetwork) <- list(tf.names, tf.names)
  if (!is.null(ppi)) {
    Idx1 <- match(ppi[, 1], tf.names)
    Idx2 <- match(ppi[, 2], tf.names)
    ok <- !is.na(Idx1) & !is.na(Idx2)
    w <- as.numeric(ppi[ok, 3])
    tfCoopNetwork[cbind(Idx1[ok], Idx2[ok])] <- w
    tfCoopNetwork[cbind(Idx2[ok], Idx1[ok])] <- w
  }
  list(expr = expr, regulatoryNetwork = regulatoryNetwork, tfCoopNetwork = tfCoopNetwork,
       tf.names = tf.names, gene.names = gene.names)
}
