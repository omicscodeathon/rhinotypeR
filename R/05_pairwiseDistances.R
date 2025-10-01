#' Compute pairwise genetic distances between aligned DNA sequences
#'
#' @description
#' Calculates pairwise genetic distances between all sequences in an aligned
#' [Biostrings::DNAStringSet-class] object using [MSA2dist::dnastring2dist()].
#' Supports global gap deletion or default pairwise gap masking.
#'
#' @param fastaData A [Biostrings::DNAStringSet-class] object containing aligned
#'   DNA sequences. Sequences must all be the same width (e.g., from [msa::msa()]
#'   or [alignToRefs()]).
#' @param model A character string specifying the substitution model to use.
#'   Defaults to `"IUPAC"`, which accounts for ambiguous bases. Other models
#'   are passed to [ape::dist.dna()] through [MSA2dist::dnastring2dist()].
#' @param deleteGapsGlobally Logical. If `TRUE`, all columns containing at least
#'   one gap (`-`) are removed prior to distance calculation (global gap deletion).
#'   If `FALSE` (default), gaps are handled pairwise during distance calculation.
#' @param ... Additional arguments passed to [ape::dist.dna()] via
#'   [MSA2dist::dnastring2dist()].
#'
#' @details
#' - If `deleteGapsGlobally = TRUE`, all gapped sites are removed across the
#'   entire alignment using [deleteMissingDataSites()].
#' - Otherwise, gaps are masked on a pairwise basis by
#'   [MSA2dist::dnastring2dist()].
#'
#' @return A numeric matrix of pairwise genetic distances, with row and column
#'   names corresponding to the sequence names in `fastaData`.
#'
#' @seealso
#' * [countSNPs()] for converting distances to SNP counts.
#' * [deleteMissingDataSites()] for global gap trimming.
#'
#' @examples
#' if (interactive()) {
#' 
#'   aln <- DNAStringSet(c(
#'     Seq1 = "ATGC",
#'     Seq2 = "ATGT",
#'     Seq3 = "ATGA"
#'   ))
#'
#'   # Default: IUPAC model with pairwise gap masking
#'   pairwiseDistances(aln)
#'   
#'   # Jukes and Cantor model
#'   pairwiseDistances(aln, model = "JC69")
#' }
#' 
#' @importFrom MSA2dist dnastring2dist
#' @export
pairwiseDistances <- function(fastaData,
                              model = "IUPAC",
                              deleteGapsGlobally = FALSE,
                              ...) {
  # checks
  if (!inherits(fastaData, "DNAStringSet"))
    stop("fastaData must be a Biostrings::DNAStringSet")
  w <- Biostrings::width(fastaData)
  if (length(unique(w)) != 1L)
    stop("Sequences must be aligned (identical widths).")
  
  if (deleteGapsGlobally) {
    fastaData <- deleteMissingDataSites(fastaData)
  }
  
  # Compute pairwise distances
  res <- MSA2dist::dnastring2dist(
    dna   = fastaData,
    model = model,   # substitution model
    ...              # other ape::dist.dna parameters
  )
  
  # Extract the distance matrix
  distanceMatrix <- as.matrix(res$distSTRING)
  
  return(distanceMatrix)
}
