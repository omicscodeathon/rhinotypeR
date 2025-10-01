#' Plot Pairwise Distance Heatmap
#'
#' Generates a heatmap from a pairwise distance matrix, typically produced by
#' \code{\link{pairwiseDistances}}.
#'
#' @param distancesMatrix A numeric matrix of pairwise distances, usually
#'   produced by \code{\link{pairwiseDistances}} or related functions.
#'
#' @details
#' The function uses the base R \code{heatmap} function to visualize distances.
#' Row and column clustering are disabled to preserve the input ordering.
#'
#' Colors are drawn from \code{heat.colors(256)}. The scale is set to "none"
#' so the distances are displayed directly, not normalized by row or column.
#'
#' @family visualization
#'
#' @seealso \code{\link{pairwiseDistances}}, \code{\link{plotTree}}
#'
#' @examples
#' # Example using built-in dataset
#' test <- system.file("extdata", "input_aln.fasta", package = "rhinotypeR")
#' fastaD <- Biostrings::readDNAStringSet(test)
#'
#' # Compute pairwise distances
#' dist_mat <- pairwiseDistances(fastaD)
#'
#' # Plot heatmap of distances
#'  plotDistances(dist_mat)
#'  
#' @return Invisibly returns the object from \code{stats::heatmap}
#'   (a list with components such as \code{rowInd} and \code{colInd}).
#'   The primary output is the heatmap drawn to the active device.
#'    
#' @importFrom stats heatmap
#' @importFrom grDevices heat.colors
#' @export
plotDistances <- function(distancesMatrix) {
  # Convert the data to a matrix 
  distancesMatrix <- as.matrix(distancesMatrix)
  
  # Heatmap
  myplot <- stats::heatmap(
    distancesMatrix,
    Rowv = NA,
    Colv = NA,
    col = heat.colors(256),
    scale = "none",
    margins = c(7, 10)
  )
  
  invisible(myplot)
}
