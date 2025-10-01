#' Overall Mean Pairwise Distance
#'
#' Computes the overall mean pairwise distance between all sequences
#' in an alignment. Optionally deletes columns containing gaps across all sequences
#' before computing the distances.
#'
#' @param fastaData A \link[Biostrings]{DNAStringSet} containing aligned nucleotide sequences.
#' @param model Character string specifying the substitution model to use.
#'   Defaults to `"IUPAC"`. Other options are passed to \code{\link[ape]{dist.dna}}.
#' @param deleteGapsGlobally Logical; if \code{TRUE}, columns with gaps in any sequence
#'   are removed before computing distances. Default is \code{FALSE}.
#' @param ... Additional arguments passed to \code{\link[MSA2dist]{dnastring2dist}}.
#'
#' @return A single numeric value giving the overall mean pairwise distance
#'   across all sequence comparisons (excluding diagonal/self comparisons).
#'
#' @details
#' The function computes all pairwise distances using
#' \code{\link[MSA2dist]{dnastring2dist}}, extracts the distance matrix,
#' and reports the mean of the off-diagonal elements.
#'
#' @seealso \code{\link{pairwiseDistances}}, \code{\link{assignTypes}}
#'
#' @examples
#' # Load example alignment
#' test <- system.file("extdata", "input_aln.fasta", package = "rhinotypeR")
#' fastaD <- Biostrings::readDNAStringSet(test)
#'
#' # Compute mean distance (using default IUPAC model)
#' overallMeanDistance(fastaD)
#'
#' # Compute mean distance with gaps deleted
#' overallMeanDistance(fastaD, deleteGapsGlobally = TRUE)
#'
#' @export
overallMeanDistance <- function(fastaData, model = "IUPAC", deleteGapsGlobally = FALSE, ...) {

  
  if (deleteGapsGlobally) {
    fastaData <- deleteMissingDataSites(fastaData)
  }
  
  # Compute pairwise distances
  res <- MSA2dist::dnastring2dist(
    dna   = fastaData,
    model = model,   # any model from ape::dist.dna [default: IUPAC]
    ...              # additional ape::dist.dna parameters
  )
  
  # Extract the distance matrix
  distanceMatrix <- as.matrix(res$distSTRING)
  
  # Overall mean pairwise distance (exclude diagonal)
  meanDist <- mean(distanceMatrix[lower.tri(distanceMatrix)])
  
  return(meanDist)
}
