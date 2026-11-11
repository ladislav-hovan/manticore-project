#' @importFrom Matrix Matrix
#' @export lioness_single
#' @title Build gene regulatory networks for a single sample using PANDA and the LIONESS equation
#' @description
#' Constructs gene regulatory networks from single-cell/nuclei RNA-seq data
#' for a single sample by running the PANDA (Passing Attributes between
#' Networks for Data Assimilation) message-passing algorithm on the expression
#' matrix without that sample to integrate transcription factor motifs,
#' protein-protein interactions, and gene expression data into unified
#' regulatory networks. The LIONESS equation is then applied to infer the
#' single-sample network for the sample of interest.
#' @author Ladislav Hovan <ladislav.hovan@ncmbm.uio.no>, Daniel Osorio <daniecos@uio.no>
#' @param pandaNet A list of networks inferred by PANDA on the full dataset (before single-sample inference).
#' @param sample The name of the sample for which to infer a single-sample network.
#' @param gexAdjusted An expression dataset, with genes in the rows and barcodes (cells) in the columns.
#' @param tfMotifs A motif dataset (data.frame or matrix) with 3 columns: TF, target gene, and motif score. Pass NULL for co-expression analysis only.
#' @param ppiNet A Protein-Protein-Interaction dataset (data.frame or matrix) with 3 columns: protein 1, protein 2, and interaction score. Pass NULL to disable protein interaction integration.
#' @param computingEngine Character specifying computing device: 'cpu' or 'gpu' (if available). Default 'cpu'.
#' @param nCores Number of processors to be used if BLAS or MPI is active.
#' @param alphaValue Numeric update parameter (0 to 1) controlling relative contribution of prior networks. Default 0.1.
#' @param hammingValue Numeric convergence threshold based on Hamming distance. Algorithm stops when updates fall below this. Default 0.001.
#' @param assocMethod Association method. Must be one of 'pearson', 'spearman' or 'pcNet'.
#' @param nIter Sets the maximum number of iterations PANDA can run before exiting.
#' @param outNet A vector containing which networks to return. Options include "regNet", "coregNet", "coopNet".
#' @param zScaling Boolean to indicate use of Z-Scores in output. FALSE will use [0,1] scale.
#' @param showProgress Boolean to indicate printing of output for algorithm progress.
#' @param randomizationMethod Method by which to randomize gene expression matrix. Default "None". Must be one of "None", "within.gene", "by.genes". "within.gene" randomization scrambles each row of the gene expression matrix, "by.gene" scrambles gene labels.
#' @param scaleByPresent Boolean to indicate scaling of correlations by percentage of positive samples.
#'   before network inference. Default FALSE.
#' @return A list of 6 elements describing the inferred networks at convergence:
#'   \itemize{
#'     \item{regNet: Regulatory network matrix (TFs × genes)}
#'     \item{coregNet: Co-regulation network matrix (genes × genes)}
#'     \item{coopNet: Cooperation network matrix (TFs × TFs)}
#'     \item{numGenes: Number of genes in the network}
#'     \item{numTFs: Number of transcription factors}
#'     \item{numEdges: Total number of edges in regulatory network}
#'   }
#' @seealso \code{\link{runMANTICORE}} for building supercell-specific networks, which uses this function internally.
lioness_single <- function(
  pandaNet,
  sample,
  gexAdjusted,
  tfMotifs = NULL,
  ppiNet = NULL,
  computingEngine = "cpu",
  nCores = 1,
  assocMethod = "pearson",
  alphaValue = 0.1,
  hammingValue = 0.001,
  nIter = Inf,
  outNet = c("regNet", "coregNet", "coopNet"),
  zScaling = TRUE,
  showProgress = TRUE,
  randomizationMethod = "None",
  scaleByPresent = FALSE
) {
  # Input validation
  if (!is.numeric(alphaValue) || alphaValue < 0 || alphaValue > 1) {
    cli::cli_abort("alphaValue must be a numeric value between 0 and 1")
  }
  if (!is.numeric(hammingValue) || hammingValue < 0) {
    cli::cli_abort("hammingValue must be a non-negative numeric value")
  }
  if (!assocMethod %in% c("pearson", "spearman", "pcNet")) {
    cli::cli_abort("assocMethod must be one of: 'pearson', 'spearman', 'pcNet'")
  }
  if (!computingEngine %in% c("cpu", "gpu")) {
    cli::cli_abort("computingEngine must be either 'cpu' or 'gpu'")
  }

  if (showProgress) {
    cli::cli_h1("MANTICORE")
  }

  # Control BLAS threading to respect nCores
  if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
    old_blas <- RhpcBLASctl::blas_get_num_procs()
    RhpcBLASctl::blas_set_num_threads(nCores)
    on.exit(RhpcBLASctl::blas_set_num_threads(old_blas), add = TRUE)
  }

  if (is.null(ppiNet) & is.null(tfMotifs)) {
    # Densify once; centering already destroys sparsity, so keep a single dense
    # base copy and use BLAS tcrossprod (dsyrk) for the co-expression matrix.
    gexAdjusted <- as.matrix(gexAdjusted)
    if (assocMethod == "spearman") {
      gexAdjusted <- t(apply(gexAdjusted, 1, rank))
    }
    geneCoExpr <- gexAdjusted - rowMeans(gexAdjusted)
    geneCoExpr <- geneCoExpr / sqrt(rowSums(geneCoExpr^2))
    return(tcrossprod(geneCoExpr))
  }
  outNetworks <- runPANDA(
    motif = tfMotifs,
    ppi = ppiNet,
    expr = gexAdjusted[, colnames(gexAdjusted) != sample, drop = FALSE],
    computing.engine = computingEngine,
    n.cores = nCores,
    alpha = alphaValue,
    hamming = hammingValue,
    iter = nIter,
    output = outNet,
    zScale = zScaling,
    progress = showProgress,
    randomize = randomizationMethod,
    assoc.method = assocMethod,
    scale.by.present = scaleByPresent
  )
  # Apply the LIONESS equation to compute single-sample network
  for (netName in names(outNetworks)) {
    if (is.matrix(outNetworks[[netName]])) {
      outNetworks[[netName]] <- outNetworks[[netName]] + (
        ncol(gexAdjusted) * (pandaNet[[netName]] - outNetworks[[netName]]))
    }
  }

  return(outNetworks)
}