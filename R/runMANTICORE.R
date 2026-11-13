#' @title Run LIONESS across cell groups and return combined networks
#' @description Builds per-supercell regulatory networks by running LIONESS on supercells from subsets of cells defined by \code{cellsMetadata} and combining the resulting networks into a wide-format data frame where each column corresponds to a network.
#' @author Ladislav Hovan <ladislav.hovan@ncmbm.uio.no>, Daniel Osorio <daniecos@uio.no>
#' @param gexMatrix An expression dataset with genes in the rows and barcodes (cells) in the columns.
#' @param tfMotifs A motif dataset, a data.frame or a matrix containing 3 columns. Each row describes a motif associated with a transcription factor (column 1) a gene (column 2) and a score (column 3).
#' @param ppiNet A Protein-Protein-Interaction dataset, a data.frame or matrix containing 3 columns. Each row describes a protein-protein interaction between transcription factor 1 (column 1), transcription factor 2 (column 2) and a score (column 3).
#' @param cellsMetadata A data.frame with cell-level metadata; must contain columns specified in \code{groupBy}.
#' @param groupBy Character vector of one or more column names in \code{cellsMetadata} to use for grouping cells into networks.
#' @param normalizeData Boolean to indicate normalization of expression data. Default TRUE performs log normalization.
#' @param removeBatchEffect Boolean to indicate batch effect correction. Default FALSE.
#' @param batch Factor or vector giving batch assignment for each cell; required if \code{removeBatchEffect = TRUE}.
#' @param minCells Minimum number of cells per group required to build a network. Default is 30.
#' @param computingEngine Either 'cpu' or 'gpu'. Passed to \code{\link{lioness_single}}.
#' @param nCores Number of processors to be used if BLAS or MPI is active.
#' @param gammaValue Graining level of data (proportion of number of single cells to super-cells). Default 10.
#' @param nPC Number of principal components to use for kNN network construction. Default 25.
#' @param assocMethod Association method. Must be one of 'pearson', 'spearman' or 'pcNet'. Default 'pearson'.
#' @param alphaValue Value to be used for update variable in PANDA. Default 0.1.
#' @param hammingValue Value at which to terminate the process based on Hamming distance. Default 0.001.
#' @param nIter Sets the maximum number of iterations PANDA can run before exiting. Default Inf.
#' @param outNet Character vector specifying which network(s) to extract. Options include "regNet", "coregNet", "coopNet". Default "regNet". When more than one network is requested, an \code{edge_type} column ("tf-target", "gene-gene", "tf-tf") is added and the networks are stacked in long format.
#' @param zScaling Boolean to indicate use of Z-Scores in output. FALSE will use [0,1] scale. Default TRUE.
#' @param showProgress Boolean to indicate printing of output for algorithm progress. Default TRUE.
#' @param randomizationMethod Method by which to randomize gene expression matrix. Default "None". Must be one of "None", "within.gene", "by.gene".
#' @param scaleByPresent Boolean to indicate scaling of correlations by percentage of positive samples. Default FALSE.
#' @param filterExpr Boolean to indicate whether or not to remove genes with 0 expression across all cells. Default FALSE.
#' @return A data.frame in wide format where rows represent TF-target pairs (union across all networks) and columns represent network identifiers. Cell values are edge weights from the corresponding network. When multiple network types are requested via \code{outNet}, an additional leading \code{edge_type} column identifies the network each row comes from and the network types are stacked in long format.
#' @seealso \code{\link{testEdges}}, \code{\link{regressEdges}}
#' @examples
#' \dontrun{
#' # Load test data
#' data(manticoreTest)
#'
#' # Example 1: Group by single column (region)
#' nets_by_region <- runMANTICORE(
#'   gexMatrix = manticoreTest$gex,
#'   tfMotifs = manticoreTest$tf,
#'   ppiNet = manticoreTest$ppi,
#'   cellsMetadata = manticoreTest$metadata,
#'   groupBy = "region"
#' )
#'
#' # -- MANTICORE -------------------------------------------------------------
#' # + Normalizing data (log scale)
#' # i 3 groups requested
#' # + 3 groups meet the minimum cell requirement (10)
#' # i Computing the background PANDA network
#' # + Initializing and validating
#' # + Verified sufficient samples
#' # i Normalizing networks
#' # i Learning network
#' # i Using tanimoto similarity
#' # i Successfully ran PANDA on 281 Genes and 963 TFs
#' # i Time elapsed: 3.99 seconds
#' # i Computing 18 networks
#' # i Networks successfully constructed
#' # i Total time elapsed: 68.53 seconds
#' # i Networks successfully combined
#'
#' # head(nets_by_region)
#' #                           tf target        T_1        T_2        T_3        T_4
#' # 1                       AATF  ACKR1 -0.2941782 -0.3354834 -0.3320515 -0.3251915
#' # 2                       ABL1  ACKR1 -0.3179234 -0.3571915 -0.3548165 -0.3481596
#' # 3                      ACSS2  ACKR1 -0.2947110 -0.3352503 -0.3322968 -0.3254164
#' # 4                       ADNP  ACKR1  0.2259662  0.2433903  0.2599372  0.2385351
#' # 5                      AEBP2  ACKR1 -0.1653191 -0.2347701 -0.1767304 -0.1722081
#' # 6 AEBP2_EED_EZH2_RBBP4_SUZ12  ACKR1 -0.2957854 -0.3360431 -0.3333447 -0.3264213
#'
#' # Example 2: Group by single column (donor)
#' nets_by_donor <- runMANTICORE(
#'   gexMatrix = manticoreTest$gex,
#'   tfMotifs = manticoreTest$tf,
#'   ppiNet = manticoreTest$ppi,
#'   cellsMetadata = manticoreTest$metadata,
#'   groupBy = "donor"
#' )
#'
#' # -- MANTICORE -------------------------------------------------------------
#' # + Normalizing data (log scale)
#' # i 3 groups requested
#' # + 3 groups meet the minimum cell requirement (10)
#' # i Computing the background PANDA network
#' # + Initializing and validating
#' # + Verified sufficient samples
#' # i Normalizing networks
#' # i Learning network
#' # i Using tanimoto similarity
#' # + Successfully ran PANDA on 281 Genes and 963 TFs
#' # i Time elapsed: 6.17 seconds
#' # i Computing 18 networks
#' # + Networks successfully constructed
#' # i Total time elapsed: 71.37 seconds
#' # + Networks successfully combined
#'
#' # head(nets_by_donor)
#' #                           tf target      P31_1      P31_2      P31_3      P31_4
#' # 1                       AATF  ACKR1 -0.3202550 -0.3172555 -0.2880045 -0.3367420
#' # 2                       ABL1  ACKR1 -0.3472954 -0.3409893 -0.3149934 -0.3576329
#' # 3                      ACSS2  ACKR1 -0.3206402 -0.3174249 -0.2885667 -0.3367325
#' # 4                       ADNP  ACKR1  0.2145595  0.2399915  0.2239300  0.2414298
#' # 5                      AEBP2  ACKR1 -0.1670806 -0.1670767 -0.1642898 -0.1716003
#' # 6 AEBP2_EED_EZH2_RBBP4_SUZ12  ACKR1 -0.3222239 -0.3183888 -0.2899739 -0.3373061
#'
#' # Example 3: Group by two columns (donor and region)
#' nets_by_donor_region <- runMANTICORE(
#'   gexMatrix = manticoreTest$gex,
#'   tfMotifs = manticoreTest$tf,
#'   ppiNet = manticoreTest$ppi,
#'   cellsMetadata = manticoreTest$metadata,
#'   groupBy = c("donor", "region")
#' )
#'
#' # -- MANTICORE -------------------------------------------------------------
#' # + Normalizing data (log scale)
#' # i 9 groups requested
#' # + 9 groups meet the minimum cell requirement (10)
#' # i Computing the background PANDA network
#' # + Initializing and validating
#' # + Verified sufficient samples
#' # i Normalizing networks
#' # i Learning network
#' # i Using tanimoto similarity
#' # + Successfully ran PANDA on 281 Genes and 963 TFs
#' # i Time elapsed: 3.77 seconds
#' # i Computing 18 networks
#' # + Networks successfully constructed
#' # i Total time elapsed: 70.28 seconds
#' # + Networks successfully combined
#' 
#' # head(nets_by_donor_region)
#' #                           tf target   P31--T_1   P31--T_2   P31--B_1   P31--B_2
#' # 1                       AATF  ACKR1 -0.3190902 -0.2185745 -0.3152367 -0.3230080
#' # 2                       ABL1  ACKR1 -0.3408164 -0.2572263 -0.3390597 -0.3470957
#' # 3                      ACSS2  ACKR1 -0.3194359 -0.2202008 -0.3154443 -0.3227614
#' # 4                       ADNP  ACKR1  0.2474420  0.1767450  0.2361739  0.2282932
#' # 5                      AEBP2  ACKR1 -0.1641631 -0.1094063 -0.1629011 -0.1673652
#' # 6 AEBP2_EED_EZH2_RBBP4_SUZ12  ACKR1 -0.3200036 -0.2212101 -0.3168823 -0.3244148
#'
#' # Example 4: Group by three columns (donor, region, and cell_type)
#' nets_by_donor_region_cell_type <- runMANTICORE(
#'   gexMatrix = manticoreTest$gex,
#'   tfMotifs = manticoreTest$tf,
#'   ppiNet = manticoreTest$ppi,
#'   cellsMetadata = manticoreTest$metadata,
#'   groupBy = c("donor", "region", "cell_type")
#' )
#'
#' # -- MANTICORE -----------------------------------
#' # + Normalizing data (log scale)
#' # i 9 groups requested
#' # + 9 groups meet the minimum cell requirement (10)
#' # i Computing the background PANDA network
#' # + Initializing and validating
#' # + Verified sufficient samples
#' # i Normalizing networks
#' # i Learning network
#' # i Using tanimoto similarity
#' # + Successfully ran PANDA on 281 Genes and 963 TFs
#' # i Time elapsed: 3.84 seconds
#' # i Computing 18 networks
#' # + Networks successfully constructed
#' # i Total time elapsed: 67.65 seconds
#' # + Networks successfully combined
#'
#' # head(nets_by_donor_region_cell_type)
#' #                           tf target P31--T--Epithelial_1 P31--T--Epithelial_2
#' # 1                       AATF  ACKR1           -0.3190902           -0.2185745
#' # 2                       ABL1  ACKR1           -0.3408164           -0.2572263
#' # 3                      ACSS2  ACKR1           -0.3194359           -0.2202008
#' # 4                       ADNP  ACKR1            0.2474420            0.1767450
#' # 5                      AEBP2  ACKR1           -0.1641631           -0.1094063
#' # 6 AEBP2_EED_EZH2_RBBP4_SUZ12  ACKR1           -0.3200036           -0.2212101
#'
#' # Example 5: Using GPU computing engine (if available)
#' nets_gpu <- runMANTICORE(
#'   gexMatrix = manticoreTest$gex,
#'   tfMotifs = manticoreTest$tf,
#'   ppiNet = manticoreTest$ppi,
#'   cellsMetadata = manticoreTest$metadata,
#'   groupBy = "region",
#'   computingEngine = "gpu"
#' )
#'
#' # -- MANTICORE -------------------------------------------------------------
#' # + Normalizing data (log scale)
#' # i 3 groups requested
#' # + 3 groups meet the minimum cell requirement (10)
#' # i Computing the background PANDA network
#' # + Initializing and validating
#' # + Verified sufficient samples
#' # i Normalizing networks
#' # i Learning network
#' # i Using tanimoto similarity
#' # + Successfully ran PANDA on 281 Genes and 963 TFs
#' # i Time elapsed: 5.14 seconds
#' # i Computing 18 networks
#' # + Networks successfully constructed
#' # i Total time elapsed: 104.15 seconds
#' # + Networks successfully combined
#'
#' # head(nets_gpu)
#' #                           tf target        T_1        T_2        T_3        T_4
#' # 1                       AATF  ACKR1 -0.3149668 -0.2882942 -0.3105499 -0.2644840
#' # 2                       ABL1  ACKR1 -0.3384307 -0.3121278 -0.3344031 -0.2886737
#' # 3                      ACSS2  ACKR1 -0.3147451 -0.2885594 -0.3110744 -0.2651064
#' # 4                       ADNP  ACKR1  0.2332610  0.2142980  0.2519617  0.2150365
#' # 5                      AEBP2  ACKR1 -0.1630593 -0.1639258 -0.1621088 -0.1211030
#' # 6 AEBP2_EED_EZH2_RBBP4_SUZ12  ACKR1 -0.3163287 -0.2902277 -0.3113921 -0.2666489
#'
#' # Example 6: Removing batch effect using donor as batch
#' nets_batch_corrected <- runMANTICORE(
#'   gexMatrix = manticoreTest$gex,
#'   tfMotifs = manticoreTest$tf,
#'   ppiNet = manticoreTest$ppi,
#'   cellsMetadata = manticoreTest$metadata,
#'   groupBy = "region",
#'   removeBatchEffect = TRUE,
#'   batch = manticoreTest$metadata$donor
#' )
#'
#' # -- MANTICORE -------------------------------------------------------------
#' # + Normalizing data (log scale)
#' # + Correcting for batch effects
#' # i 3 groups requested
#' # + 3 groups meet the minimum cell requirement (10)
#' # i Computing the background PANDA network
#' # + Initializing and validating
#' # + Verified sufficient samples
#' # i Normalizing networks
#' # i Learning network
#' # i Using tanimoto similarity
#' # + Successfully ran PANDA on 281 Genes and 963 TFs
#' # i Time elapsed: 3.66 seconds
#' # i Computing 18 networks
#' # + Networks successfully constructed
#' # i Total time elapsed: 66.84 seconds
#' # + Networks successfully combined
#'
#' # head(nets_batch_corrected)
#' #                           tf target        T_1        T_2        T_3        T_4
#' # 1                       AATF  ACKR1 -0.2962846 -0.3274905 -0.2734708 -0.3312162
#' # 2                       ABL1  ACKR1 -0.3224845 -0.3490377 -0.3011509 -0.3564239
#' # 3                      ACSS2  ACKR1 -0.2971183 -0.3278750 -0.2746225 -0.3306292
#' # 4                       ADNP  ACKR1  0.2216245  0.2543183  0.2191569  0.2430913
#' # 5                      AEBP2  ACKR1 -0.1643446 -0.0555660 -0.1490343 -0.1655627
#' # 6 AEBP2_EED_EZH2_RBBP4_SUZ12  ACKR1 -0.2991352 -0.3258038 -0.2779301 -0.3312577
#' }
#' @export
#' @importFrom dplyr %>% mutate group_by filter .data pull summarize
#' @importFrom cli cli_h1 cli_alert_success cli_alert_info cli_abort cli_progress_along
#' @importFrom furrr future_map furrr_options
#' @importFrom future plan multisession sequential
runMANTICORE <- function(
  gexMatrix,
  tfMotifs,
  ppiNet,
  cellsMetadata,
  groupBy,
  normalizeData = TRUE,
  removeBatchEffect = FALSE,
  batch = NULL,
  minCells = 30,
  computingEngine = "cpu",
  nCores = 1,
  gammaValue = 10,
  nPC = 25,
  assocMethod = "pearson",
  alphaValue = 0.1,
  hammingValue = 0.001,
  nIter = Inf,
  outNet = "regNet",
  zScaling = TRUE,
  showProgress = TRUE,
  randomizationMethod = "None",
  scaleByPresent = FALSE,
  filterExpr = FALSE
) {
  # Input validation
  if (ncol(gexMatrix) != nrow(cellsMetadata)) {
    cli::cli_abort("gexMatrix must have the same number of columns as cellsMetadata has rows")
  }
  if (!all(groupBy %in% colnames(cellsMetadata))) {
    cli::cli_abort("groupBy columns not found in cellsMetadata: {paste(setdiff(groupBy, colnames(cellsMetadata)), collapse=', ')}")
  }
  if (!is.numeric(minCells) || minCells < 1) {
    cli::cli_abort("minCells must be a positive integer")
  }
  if (removeBatchEffect && is.null(batch)) {
    cli::cli_abort("batch must be provided when removeBatchEffect = TRUE")
  }
  if (!is.null(batch) && length(batch) != ncol(gexMatrix)) {
    cli::cli_abort("batch must have the same length as gexMatrix columns (cells)")
  }
  if (alphaValue < 0 || alphaValue > 1) {
    cli::cli_abort("alphaValue must be a numeric value between 0 and 1")
  }
  validNets <- c("regNet", "coregNet", "coopNet")
  if (length(outNet) < 1 || !all(outNet %in% validNets)) {
    cli::cli_abort("outNet must be one or more of: {paste(validNets, collapse=', ')}")
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

  # Normalizing data
  if (normalizeData) {
    if (showProgress) {
      cli::cli_alert_success("Normalizing data (log scale)")
    }
    gexMatrix <- log_normalize_data(gexMatrix)
  }

  # Removing batch effect
  if (removeBatchEffect) {
    if (is.null(batch)) {
      cli::cli_abort('batch is needed for batch effect correction')
    }
    if (length(unique(batch[!is.na(batch)])) < 2) {
      # A single batch leaves nothing to correct; skip so the per-gene median
      # is not added to uncorrected data.
      cli::cli_alert_warning("batch has fewer than two levels; skipping batch effect correction")
    } else {
      if (showProgress) {
        cli::cli_alert_success("Correcting for batch effects")
      }
      mean_expr <- apply(gexMatrix, 1, median)
      gexMatrix <- remove_batch(X = gexMatrix, batch = batch)
      gexMatrix <- gexMatrix + mean_expr
    }
  }
  rm(batch)

  # Pre-convert to data.frame once so workers receive the converted objects
  # instead of converting on every call
  tfMotifs <- as.data.frame(tfMotifs)
  ppiNet <- as.data.frame(ppiNet)

  # Setting min number of cells to construct network
  min_cells <- max(minCells, 3 * gammaValue)  # Ensure at least 3 super-cells
  min_cells_per_group <- gammaValue

  if (all(groupBy %in% colnames(cellsMetadata))) {
    # Vectorized row-wise concatenation of the grouping columns; avoids the
    # per-row apply() closure over cellsMetadata.
    collapsedGroup <- do.call(paste, c(cellsMetadata[, groupBy, drop = FALSE], sep = "--"))
    metadata <- data.frame(cell_id = colnames(gexMatrix), network_id = collapsedGroup)
    metadata <- metadata %>%
      group_by(.data$network_id) %>%
      mutate(n_cells = length(.data$cell_id))
  } else {
    cli::cli_abort('groupBy must match cellsMetadata column name')
  }
  rm(cellsMetadata)

  total_net <- length(unique(metadata$network_id))
  # Every subgroup needs at least 10 (by default) cells
  metadata <- metadata %>% filter(.data$n_cells >= min_cells_per_group)
  filtered_net <- length(unique(metadata$network_id))

  if (showProgress) {
    cli::cli_alert_info(paste0(total_net, " groups requested"))
    cli::cli_alert_success(paste0(filtered_net, " groups meet the minimum cell requirement (", min_cells_per_group, ")"))
  }

  # Check if there are enough total cells remaining after filtering to build reliable networks (by default 30)
  n_total_cells <- sum(metadata %>% group_by(.data$network_id) %>% summarize(mean_n_cells = mean(.data$n_cells)) %>% pull(.data$mean_n_cells))
  if (n_total_cells < min_cells) {
    cli::cli_abort(paste0("There are not enough total cells (< ", min_cells, ") to build reliable networks."))
  }

  network_ids <- unique(metadata$network_id)
  cells_by_network <- split(metadata$cell_id, metadata$network_id)

  # Group the supercells here but keep everything together
  gexList = list()
  for (nid in network_ids) {
    gexChunk <- gexMatrix[, cells_by_network[[nid]], drop = FALSE]
    gexSC <- makeSuperCells(X = gexChunk, gamma = gammaValue, n.pc = nPC, fast.pca = FALSE)
    n_cols <- ncol(gexSC)
    colnames(gexSC) <- paste0(nid, "_", seq_len(n_cols))
    gexList[[nid]] <- gexSC
  }
  gexAdjusted <- do.call(cbind, gexList)
  rm(gexMatrix, metadata, gexList, gexChunk, gexSC, cells_by_network)

  if (showProgress) {
    cli::cli_alert_info(paste0("Computing the background PANDA network"))
  }
  pandaNet <- runPANDA(
    motif = tfMotifs,
    ppi = ppiNet,
    expr = gexAdjusted,
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

  # Applies LIONESS to an expression chunk with a single missing column
  compute_network <- function(sample) {
    network <- lioness_single(
      gexAdjusted = gexAdjusted,
      tfMotifs = tfMotifs,
      ppiNet = ppiNet,
      pandaNet = pandaNet,
      sample = sample,
      computingEngine = computingEngine,
      nCores = 1,
      assocMethod = assocMethod,
      alphaValue = alphaValue,
      hammingValue = hammingValue,
      nIter = nIter,
      outNet = outNet,
      zScaling = zScaling,
      showProgress = FALSE,
      randomizationMethod = randomizationMethod,
      scaleByPresent = scaleByPresent
    )

    # Return a single matrix when one network is requested, otherwise a
    # named list of the requested network matrices.
    if (length(outNet) == 1L) {
      return(network[[outNet]])
    }
    return(network[outNet])
  }

  # The TF-motif prior and PPI network are broadcast to every parallel worker.
  # If the optional 'mori' package is available, place them in OS-backed shared
  # memory so workers map the same physical pages instead of each receiving a
  # full serialized copy. Falls back to standard serialization when absent, or
  # when disabled via options(scorpion.use_mori = FALSE).
  use_mori <- nCores > 1 &&
    isTRUE(getOption("scorpion.use_mori", TRUE)) &&
    requireNamespace("mori", quietly = TRUE)
  if (use_mori) {
    if (!is.null(tfMotifs)) tfMotifs <- mori::share(tfMotifs)
    if (!is.null(ppiNet)) ppiNet <- mori::share(ppiNet)
    gexAdjusted <- mori::share(gexAdjusted)
    pandaNet <- mori::share(pandaNet)
  }

  furrr_opts <- furrr::furrr_options(
    seed = TRUE,
    packages = if (use_mori) "mori" else NULL
  )

  if (nCores > 1) {
    # Allow arbitrarily large globals to be exported to workers.
    old_maxsize <- getOption("future.globals.maxSize")
    options(future.globals.maxSize = Inf)
    on.exit(options(future.globals.maxSize = old_maxsize), add = TRUE)
    old_plan <- future::plan(future::multisession, workers = nCores)
    on.exit(future::plan(old_plan), add = TRUE)
  } else {
    old_plan <- future::plan(future::sequential)
    on.exit(future::plan(old_plan), add = TRUE)
  }

  n_total <- ncol(gexAdjusted)
  col_names <- colnames(gexAdjusted)

  if (showProgress) {
    cli::cli_alert_info(paste0("Computing ", n_total, " networks"))
    tic <- proc.time()[3]
    if (nCores > 1) {
      cli::cli_alert_info(paste0("Using ", nCores, " cores for parallel processing"))
      network_matrices <- furrr::future_map(col_names, compute_network, .options = furrr_opts, .progress = FALSE)
    } else {
      network_matrices <- vector("list", n_total)
      for (i in seq_len(n_total)) {
        cli::cli_alert_info(paste0("Network ", i, "/", n_total, ": ", col_names[i]))
        network_matrices[[i]] <- compute_network(col_names[i])
      }
    }
    cli::cli_alert_success("Networks successfully constructed")
    toc <- proc.time()[3]
    cli::cli_alert_info(paste0("Total time elapsed: ", round(toc - tic, 2), " seconds"))
  } else {
    network_matrices <- furrr::future_map(col_names, compute_network, .options = furrr_opts, .progress = FALSE)
  }

  n_nets <- length(network_matrices)

  if (length(outNet) == 1L) {
    # Single network type: wide format with tf, target and one column per group.
    # Per-group networks can differ in gene set/ordering (e.g. when
    # filterExpr = TRUE drops different zero-expression genes per group), so
    # build the tf/target grid from the union of all networks' dimnames
    # instead of assuming they all match the first network's shape.
    all_tfs <- unique(unlist(lapply(network_matrices, rownames)))
    all_targets <- unique(unlist(lapply(network_matrices, colnames)))
    n_tf <- length(all_tfs)
    n_target <- length(all_targets)
    tf_target_df <- expand.grid(
      tf = all_tfs,
      target = all_targets,
      KEEP.OUT.ATTRS = FALSE,
      stringsAsFactors = FALSE
    )
    n_edges <- n_tf * n_target
    # Stream extraction: place each network's weights into their global
    # tf/target positions, then NULL out the list element immediately to
    # free memory.
    weight_matrix <- matrix(NA_real_, nrow = n_edges, ncol = n_nets)
    colnames(weight_matrix) <- col_names
    for (k in seq_len(n_nets)) {
      net_k <- as.matrix(network_matrices[[k]])
      row_idx <- match(rownames(net_k), all_tfs)
      col_idx <- match(colnames(net_k), all_targets)
      global_idx <- outer(row_idx, (col_idx - 1L) * n_tf, FUN = "+")
      weight_matrix[as.vector(global_idx), k] <- as.vector(net_k)
      network_matrices[k] <- list(NULL)
    }
    rm(network_matrices)
    # Combine into final data frame
    networks <- data.frame(
      tf = as.character(tf_target_df$tf),
      target = as.character(tf_target_df$target),
      weight_matrix,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  } else {
    # Multiple network types: stack each type in long format and prepend an
    # edge_type column identifying the network each block of rows comes from.
    edge_type_labels <- c(
      regNet = "tf-target",
      coregNet = "gene-gene",
      coopNet = "tf-tf"
    )

    per_type <- lapply(outNet, function(net_name) {
      # Node pairs are specific to each network type (TFs x genes, genes x
      # genes or TFs x TFs). Per-group networks can differ in gene set/order
      # (e.g. filterExpr = TRUE drops different zero-expression genes per
      # group), so build the grid from the union of all groups' dimnames.
      all_tfs <- unique(unlist(lapply(network_matrices, function(net) rownames(net[[net_name]]))))
      all_targets <- unique(unlist(lapply(network_matrices, function(net) colnames(net[[net_name]]))))
      n_tf <- length(all_tfs)
      n_target <- length(all_targets)
      pair_df <- expand.grid(
        tf = all_tfs,
        target = all_targets,
        KEEP.OUT.ATTRS = FALSE,
        stringsAsFactors = FALSE
      )
      n_edges <- n_tf * n_target

      weight_matrix <- matrix(NA_real_, nrow = n_edges, ncol = n_nets)
      colnames(weight_matrix) <- network_ids
      for (k in seq_len(n_nets)) {
        net_k <- as.matrix(network_matrices[[k]][[net_name]])
        row_idx <- match(rownames(net_k), all_tfs)
        col_idx <- match(colnames(net_k), all_targets)
        global_idx <- outer(row_idx, (col_idx - 1L) * n_tf, FUN = "+")
        weight_matrix[as.vector(global_idx), k] <- as.vector(net_k)
      }

      data.frame(
        edge_type = edge_type_labels[[net_name]],
        tf = as.character(pair_df$tf),
        target = as.character(pair_df$target),
        weight_matrix,
        stringsAsFactors = FALSE,
        check.names = FALSE
      )
    })
    rm(network_matrices)

    networks <- do.call(rbind, per_type)
    rownames(networks) <- NULL
  }

  if (showProgress) {
    cli::cli_alert_success("Networks successfully combined")
  }

  return(networks)
}