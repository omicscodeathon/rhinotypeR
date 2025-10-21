#' Align user sequences to packaged rhinovirus prototype references
#'
#' Combines imported sequences with the prototype rhinovirus references shipped 
#' in \pkg{rhinotypeR}. Performs a multiple sequence alignment using ClustalW, 
#' ClustalOmega, or MUSCLE. The aligned result is returned as a 
#' \code{DNAStringSet} with equal-width sequences, suitable for downstream 
#' distance calculations (e.g., \code{\link{pairwiseDistances}}). 
#' If \code{trimToRef = TRUE}, the final alignment is cropped to the 
#' non-gap span of a chosen reference sequence. 
#' 
#' @param seqData A \code{Biostrings::DNAStringSet} containing the user
#'   sequences to be aligned.
#' @param method Character; one of \code{"ClustalW"}, \code{"ClustalOmega"},
#'   or \code{"Muscle"}. Passed to \code{msa::msa}.
#' @param trimToRef \code{logical(1)}. If \code{TRUE} (default), crop the aligned
#'   \code{DNAStringSet} to the span of the reference sequence (useful when user
#'   sequences extend beyond the targeted VP4/2 region).
#' @param refName \code{character(1)}. Name of the reference sequence that
#'   defines the trimming span when \code{trimToRef = TRUE}. Defaults to
#'   \code{"JN855971.1_A107"} (a VP4/2 prototype used as a convenient anchor;
#'   it has no special biological significance). Provide a different accession
#'   to change the trimming anchor.
#' @param ... Additional arguments forwarded to \code{msa::msa}
#'   (e.g., gap penalties).
#'
#' @details
#' Packaged prototypes are loaded from
#' \code{system.file("extdata","prototypes.fasta", package = "rhinotypeR")}.
#' User sequences are appended and jointly aligned with \code{msa::msa}. The
#' returned \code{DNAStringSet} has equal width. When \code{trimToRef = TRUE},
#' the alignment is cropped to the range covered by non-gap positions in
#' \code{refName}; internal gaps within that span are preserved. For global
#' removal of gap-containing columns across all sequences, see
#' \code{deleteMissingDataSites()}.
#'
#' If you prefer to align/curate outside R, use
#' \code{\link{getPrototypeSeqs}} to copy the prototype references to disk and
#' import your curated alignment later.
#'
#' @section Notes / Caveats:
#' \itemize{
#'   \item \strong{Method selection:} \code{msa::msa} expects a single method
#'         string; this wrapper enforces that with \code{match.arg()}.
#'   \item \strong{Homology:} Ensure sequences are homologous and from the same
#'         genomic region; MSA quality degrades on mixed loci.
#'   \item \strong{Result class:} The function returns a \code{DNAStringSet}
#'         (not a \code{DNAMultipleAlignment}).
#' }
#'
#' @return A \code{Biostrings::DNAStringSet} of aligned sequences (user + prototypes),
#'   all with identical width; if \code{trimToRef = TRUE}, cropped to the chosen
#'   reference’s non-gap span.
#'
#' @seealso \code{\link{getPrototypeSeqs}}, \code{\link{pairwiseDistances}},
#'   \code{\link{assignTypes}}
#'
#' @examples
#' seqs_path <- system.file("extdata", "test.fasta", package = "rhinotypeR")
#' seqs <- Biostrings::readDNAStringSet(seqs_path)
#'
#' aln_trim <- alignToRefs(seqs, method = "ClustalW", trimToRef = TRUE)
#'
#' # Keep full joint alignment (no trimming)
#' aln_full <- alignToRefs(seqs, method = "ClustalW", trimToRef = FALSE)
#'
#' # Use a different prototype as trimming anchor
#' # aln_sel <- alignToRefs(seqs, method = "Muscle", trimToRef = TRUE,
#' #                        refName = "AF343653.1_B27")
#'
#' @importFrom Biostrings readDNAStringSet DNAStringSet subseq
#' @importFrom msa msa
#' @importFrom methods as
#' @importClassesFrom Biostrings DNAMultipleAlignment
#' @export
alignToRefs <- function(seqData,
                        method = c("ClustalW", "ClustalOmega", "Muscle"),
                        trimToRef = TRUE,
                        refName = "JN855971.1_A107",
                        ...) {
  # validate input
  method <- match.arg(method, c("ClustalW", "ClustalOmega", "Muscle"))
  
  # get references
  ref <- system.file("extdata", "prototypes.fasta", package = "rhinotypeR")
  refseqs <- Biostrings::readDNAStringSet(ref)
  
  # Concatenate with rhinovirus references within the package
  allseqs <- c(seqData, refseqs)
  
  # Align with msa
  aln <- msa::msa(allseqs, method = method, ...)
  
  # Convert to Biostrings DNAMultipleAlignment -> DNAStringSet
  aln_bio <- as(aln, "DNAMultipleAlignment")
  aln_set <- Biostrings::DNAStringSet(aln_bio)
  
  # Optionally trim to reference non-gap span (default anchor kept unless user overrides)
  if (isTRUE(trimToRef)) {
    if (!refName %in% names(aln_set)) {
      stop("refName not found in alignment: '", refName,
           "'. Available names include e.g.: ",
           paste(utils::head(names(aln_set), 5L), collapse = ", "), " ...")
    }
    aln_set <- preProcessFastaStringSet(aln_set, ref_name = refName)
  }
  
  aln_set
}

# ---- internal helper ----
#' Trim an aligned DNAStringSet to the non-gap span of a reference
#'
#' Internal helper used by [alignToRefs()] to crop an aligned
#' \code{DNAStringSet} to the region spanned by non-gap characters in a chosen
#' reference sequence. Gaps *inside* that span are retained. To remove
#' gap-containing columns across all sequences, see [deleteMissingDataSites()].
#'
#' @param aln_set A \code{Biostrings::DNAStringSet} with equal widths.
#' @param ref_name \code{character(1)}. The sequence name that defines the
#'   trimming span (must be present in \code{names(aln_set)}). Defaults to
#'   \code{"JN855971.1_A107"}.
#'
#' @return A trimmed \code{Biostrings::DNAStringSet}.
#'
#' @keywords internal
#' @noRd
preProcessFastaStringSet <- function(aln_set, ref_name = "JN855971.1_A107") {
  stopifnot(inherits(aln_set, "DNAStringSet"))
  if (!all(Biostrings::width(aln_set) == Biostrings::width(aln_set)[1]))
    stop("Input must be an aligned DNAStringSet (all sequences same width).")
  if (!ref_name %in% names(aln_set))
    stop("Reference name not found: ", ref_name)
  
  ref <- aln_set[ref_name]
  ref_chars <- strsplit(as.character(ref), "", fixed = TRUE)[[1]]
  non_gap_idx <- which(ref_chars != "-")
  if (length(non_gap_idx) == 0)
    stop("Reference contains only gaps within the alignment.")
  
  start_pos <- non_gap_idx[1]
  end_pos   <- non_gap_idx[length(non_gap_idx)]
  
  Biostrings::subseq(aln_set, start = start_pos, end = end_pos)
}
