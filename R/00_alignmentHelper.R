#' Remove gap-containing sites from an alignment (global deletion)
#'
#' Deletes every alignment column that contains at least one gap (`-`),
#' returning an ungapped alignment of equal width across sequences. This is
#' a **global** deletion strategy (contrast with *pairwise deletion*, which
#' ignores gaps only in the sequences being compared).
#'
#' @param aln A [Biostrings::DNAStringSet-class] containing an **aligned**
#'   nucleotide alignment with equal-width sequences.
#'
#' @details
#' Use global deletion when you want all downstream calculations to avoid any
#' columns with gaps across the dataset. If you prefer pairwise deletion,
#' select it in your distance routine (e.g., via [pairwiseDistances()]).
#'
#' @return A [Biostrings::DNAStringSet-class] with gap-containing columns removed.
#'
#' @seealso
#'   [pairwiseDistances()] for distance calculation (with optional global deletion),
#'   [assignTypes()] for genotype assignment,
#'   [alignToRefs()] to generate aligned inputs.
#'
#' @examples
#' aln <- Biostrings::DNAStringSet(c(
#'   Seq1 = "ATGC-",
#'   Seq2 = "AT-CG",
#'   Seq3 = "ATGCG"
#' ))
#' deleteMissingDataSites(aln)
#'
#' @importFrom Biostrings DNAStringSet width
#' @keywords internal
#' @noRd
deleteMissingDataSites <- function(aln){
  if (!inherits(aln, "DNAStringSet")) {
    stop("Input must be a Biostrings::DNAStringSet.")
  }
  w <- Biostrings::width(aln)
  if (length(unique(w)) != 1L) {
    stop("All sequences must have the same width.")
  }
  
  gap_chars <- c("-", ".") # some aligners use "." as gap
  s <- as.character(aln)
  L <- w[1]
  M <- matrix(strsplit(paste0(s, collapse=""), "")[[1]], nrow=length(s), byrow=TRUE)
  stopifnot(ncol(M) == L)
  keep_cols <- apply(M, 2, function(col) !any(col %in% gap_chars))
  trimmed <- Biostrings::DNAStringSet(apply(M[, keep_cols, drop=FALSE], 1, paste0, collapse=""))
  names(trimmed) <- names(aln)
  trimmed
  
}
