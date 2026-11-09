#' @name manticoreTest
#' @docType data
#' @title Example single-cell colorectal cancer data for MANTICORE
#' @description
#' A list bundling the inputs required to build and compare gene regulatory
#' networks with MANTICORE, derived from a colorectal cancer single-cell
#' RNA-sequencing experiment. It contains a gene expression matrix, a
#' transcription factor motif prior, a protein-protein interaction prior, and
#' cell-level metadata.
#' @usage data(manticoreTest)
#' @format A named list with four elements:
#' \describe{
#'   \item{\code{gex}}{A \code{dgCMatrix} gene expression matrix with 300 genes
#'     (rows) and 1,954 cells (columns).}
#'   \item{\code{tf}}{A \code{data.frame} of transcription factor-target motif
#'     pairs from DoRothEA with columns \code{source_genesymbol},
#'     \code{target_genesymbol} and \code{weight} (371,738 rows).}
#'   \item{\code{ppi}}{A \code{data.frame} of protein-protein interactions with
#'     columns \code{source_genesymbol}, \code{target_genesymbol} and
#'     \code{weight} (4,076 rows).}
#'   \item{\code{metadata}}{A \code{data.frame} of cell-level annotations with
#'     columns \code{cell_id}, \code{donor}, \code{region} and \code{cell_type}
#'     (1,954 rows). Region codes are \code{T} (tumor), \code{B} (border) and
#'     \code{N} (normal).}
#' }
#' @examples
#' # Loading example data
#' data(manticoreTest)
#'
#' # The structure of the data
#' str(manticoreTest)
NULL