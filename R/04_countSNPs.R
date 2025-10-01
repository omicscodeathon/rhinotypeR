#' Count pairwise SNPs between aligned DNA sequences
#'
#' @description
#' Computes pairwise single nucleotide polymorphism (SNP) counts between all
#' sequences in an aligned [Biostrings::DNAStringSet-class] object.
#' Proportions are computed with [MSA2dist::dnastring2dist()] using the `"IUPAC"`
#' model (accounts for ambiguous bases), then converted to integer SNP counts.
#'
#' @param fastaData A [Biostrings::DNAStringSet-class] of aligned DNA sequences.
#'   Sequences must have identical width.
#' @param deleteGapsGlobally Logical. If `TRUE`, sites containing gaps (`-` or `.`)
#'   are removed across all sequences before SNP counting. Default: `FALSE`.
#'
#' @return A numeric matrix of SNP counts with sequence names as dimnames.
#'
#' @seealso [pairwiseDistances()] for distances; [SNPeek()], [plotAA()] for viz.
#'
#' @examples
#' if (interactive()) {
#'   aln <- Biostrings::DNAStringSet(c(Seq1="ATGC", Seq2="ATGT", Seq3="ATGA"))
#'   countSNPs(aln)
#'   countSNPs(aln, deleteGapsGlobally = TRUE)
#' }
#'
#' @importFrom MSA2dist dnastring2dist
#' @export
countSNPs <- function(fastaData, deleteGapsGlobally = FALSE) {
  if (!inherits(fastaData, "DNAStringSet"))
    stop("fastaData must be a Biostrings::DNAStringSet")
  w <- Biostrings::width(fastaData)
  if (length(unique(w)) != 1L)
    stop("Sequences must be aligned (identical widths).")
  
  if (deleteGapsGlobally) {
    fastaData <- deleteMissingDataSites(fastaData)
  }
  
  res <- MSA2dist::dnastring2dist(fastaData, model = "IUPAC")  # proportions
  D <- as.matrix(res$distSTRING)  # mismatch proportions
  N <- as.matrix(res$sitesUsed)   # sites used per pair
  
  # Round to nearest integer with +0.5 trick (safe for non-negative values)
  SNP <- matrix(as.integer(D * N + 0.5),
                nrow = nrow(D), dimnames = dimnames(D))
  SNP
}
