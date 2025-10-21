#' Count pairwise SNPs between aligned DNA sequences
#'
#' @description
#' Computes the absolute number of single-nucleotide differences (SNPs)
#' between all pairs of sequences in an aligned
#' [Biostrings::DNAStringSet-class] object. Optionally removes any
#' alignment columns that contain gaps or ambiguous bases before counting.
#'
#' @param fastaData A [Biostrings::DNAStringSet-class] of aligned DNA sequences.
#'   Sequences must have identical width.
#' @param deleteGapsGlobally Logical. If `TRUE`, sites containing gaps (`-` or `.`)
#'   or ambiguous bases (`"N"`, `"?"`) are removed across all sequences 
#'   before SNP counting. Default: `FALSE`.
#'
#' @return A square integer matrix of SNP counts with sequence names as 
#'    row and column dimnames.
#'
#' @seealso [pairwiseDistances()] for distances; [SNPeek()], [plotAA()] for viz.
#'
#' @examples
#' if (interactive()) {
#'   aln <- Biostrings::DNAStringSet(c(Seq1="ATGC", Seq2="ATGT", Seq3="ATGA"))
#'   countSNPs(aln)
#'   countSNPs(aln, deleteGapsGlobally = FALSE)
#' }
#'
#' @importFrom Biostrings width
#' @export
countSNPs <- function(fastaData, deleteGapsGlobally = FALSE) {
  if (!inherits(fastaData, "DNAStringSet"))
    stop("Input must be a Biostrings::DNAStringSet")
  
  # Ensure sequences are aligned
  seq_widths <- Biostrings::width(fastaData)
  if (length(unique(seq_widths)) != 1L)
    stop("Sequences must be aligned (identical lengths)")
  
  # Convert each DNA sequence into a character vector of bases
  seq_matrix <- do.call(
    rbind,
    lapply(as.character(fastaData), function(x) strsplit(x, "")[[1]])
  )
  
  # Optionally remove columns containing gaps or ambiguous bases
  if (deleteGapsGlobally) {
    keep_cols <- apply(seq_matrix, 2, function(col)
      !any(col %in% c("-", "N", "?")))
    seq_matrix <- seq_matrix[, keep_cols, drop = FALSE]
  }
  
  n <- nrow(seq_matrix)
  snp_matrix <- matrix(0L, n, n, dimnames = list(names(fastaData), names(fastaData)))
  
  # Pairwise SNP counting
  for (i in seq_len(n)) {
    for (j in seq_len(i)) {
      mism <- sum(seq_matrix[i, ] != seq_matrix[j, ], na.rm = TRUE)
      snp_matrix[i, j] <- snp_matrix[j, i] <- mism
    }
  }
  
  snp_matrix
}

