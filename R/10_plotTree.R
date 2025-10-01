#' Plot a Phylogenetic Tree from a Distance Matrix
#'
#' @description
#' Generates a simple phylogenetic tree/dendrogram from a pairwise distance
#' matrix using hierarchical clustering (complete linkage by default).
#'
#' @param distance_matrix A numeric, symmetric matrix of pairwise distances
#'   (zeros on the diagonal).
#' @param ... Additional arguments passed to \code{\link[stats]{plot.hclust}}.
#'
#' @details
#' The function converts the matrix to a \code{\link[stats]{dist}} object,
#' performs \code{\link[stats]{hclust}} with \code{method = "complete"}, and
#' plots the resulting dendrogram.
#'
#' @return Invisibly returns the \code{\link[stats]{hclust}} object.
#'
#' @family visualization
#' @seealso \code{\link{pairwiseDistances}}, \code{\link{plotDistances}}
#'
#' @examples
#' # Example using built-in dataset
#' test <- system.file("extdata", "input_aln.fasta", package = "rhinotypeR")
#' fastaD <- Biostrings::readDNAStringSet(test)
#'
#' # Compute pairwise distances
#' dist_mat <- pairwiseDistances(fastaD)
#'
#' # Plot tree
#' plotTree(dist_mat, lwd = 2)
#'
#' @importFrom stats as.dist hclust
#' @export
plotTree <- function(distance_matrix, ...) {
  # coerce and validate
  distance_matrix <- as.matrix(distance_matrix)
  if (!is.numeric(distance_matrix)) {
    stop("distance_matrix must be a numeric matrix.")
  }
  if (nrow(distance_matrix) != ncol(distance_matrix)) {
    stop("distance_matrix must be square.")
  }
  # require symmetry (tolerance to small FP noise)
  if (!isTRUE(all.equal(distance_matrix, t(distance_matrix), tolerance = 1e-12))) {
    stop("distance_matrix must be symmetric.")
  }
  # diagonal should be 0
  if (any(diag(distance_matrix) != 0, na.rm = TRUE)) {
    stop("Diagonal entries of distance_matrix must be zero.")
  }
  # no NAs
  if (anyNA(distance_matrix)) {
    stop("distance_matrix contains NA values; please remove or impute them.")
  }
  
  # build and plot
  d <- stats::as.dist(distance_matrix)
  hc <- stats::hclust(d, method = "complete")
  plot(hc, ...)  # stats::plot.hclust is the S3 method dispatched here
  invisible(hc)
}
