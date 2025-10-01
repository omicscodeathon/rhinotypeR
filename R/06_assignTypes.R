##' Assign rhinovirus genotypes based on pairwise distances to prototype strains
#'
#' @description
#' Compares query rhinovirus sequences to prototype reference strains using
#' pairwise genetic distances and assigns a genotype when the nearest distance
#' is below a user-specified threshold. If the nearest distance exceeds the
#' threshold, the sequence is returned as `"unassigned"` **but the actual nearest
#' distance and prototype reference are still reported**.
#'
#' @param fastaData A [Biostrings::DNAStringSet-class] containing an **aligned**
#'   VP4/2 nucleotide alignment. The alignment must include the rhinovirus
#'   prototype sequences (same accessions as shipped with the package). If your
#'   alignment does not already include prototypes, you have two options:
#'   \itemize{
#'     \item \emph{In R (no download):} call [alignToRefs()] to align your
#'           sequences jointly with the packaged prototypes, then pass the result
#'           here.
#'     \item \emph{External tools:} use [getPrototypeSeqs()] to export the
#'           prototypes, combine with your sequences, align/curate externally,
#'           then import the curated alignment.
#'   }
#' @param model Character string specifying the substitution model used by
#'   [pairwiseDistances()] for distance calculation. Defaults to `"IUPAC"`.
#' @param deleteGapsGlobally Logical. If `TRUE`, sites containing gaps (`-`)
#'   are removed across all sequences before distance calculation. Default: `FALSE`.
#' @param threshold Numeric. Distance threshold for genotype assignment.
#'   Defaults to `0.105` (≈10.5%), consistent with species A/C thresholds in
#'   McIntyre et al. (2013). Use `0.095` for species B if desired.
#'
#' @details
#' - Prototype accessions are distributed with the package and loaded from
#'   `inst/extdata/prototypes.rda` (object `prototypes`).
#' - The function checks that all prototypes are present in `fastaData`. If not,
#'   it stops with guidance to use [alignToRefs()] or [getPrototypeSeqs()].
#' - Distances are computed via [pairwiseDistances()] on the provided alignment.
#' - Above-threshold queries are labeled `"unassigned"` while still reporting
#'   the nearest prototype and its distance.
#'
#' @return A `data.frame` with columns:
#' \describe{
#'   \item{query}{Query sequence name.}
#'   \item{assignedType}{Assigned rhinovirus type, or `"unassigned"`.}
#'   \item{distance}{Nearest pairwise distance to any prototype (reported even when unassigned).}
#'   \item{reference}{Prototype accession corresponding to the nearest distance.}
#' }
#'
#' @seealso
#' [alignToRefs()] for in-R alignment with packaged prototypes;
#' [getPrototypeSeqs()] to export prototypes for external alignment;
#' [pairwiseDistances()] for distance calculation.
#'
#' @references
#' McIntyre, C. L., et al. (2013). Proposals for the classification of human
#' rhinovirus species A, B and C into genotypically assigned types.
#' \emph{Journal of General Virology}, 94(8), 1791–1806.
#'
#' @examples
#' if (interactive()) {
#'
#'   # Load example alignment shipped with the package
#'   test <- system.file("extdata", "input_aln.fasta", package = "rhinotypeR")
#'   fastaD <- Biostrings::readDNAStringSet(test)
#'
#'   # Assign types
#'   # Note: The input must include the prototype sequences.
#'   # If your alignment does not already contain prototypes,
#'   # use alignToRefs() before calling assignTypes() OR
#'   # use getPrototypeSeqs() to download the prototype sequences, combine with your new sequences and align before importing to R
#'   try({
#'     res <- assignTypes(fastaD, model = "IUPAC", deleteGapsGlobally = FALSE, threshold = 0.105)
#'     head(res)
#'   })
#' }
#' 
#' @export
assignTypes <- function(fastaData,
                       model = "IUPAC",
                       deleteGapsGlobally = FALSE,
                       threshold = 0.105) {
  
  if (deleteGapsGlobally) {
    fastaData <- deleteMissingDataSites(fastaData)
  }
  
  # Load prototype accessions (expects object 'prototypes' with $Accession)
  ref <- system.file("extdata", "prototypes.rda", package = "rhinotypeR")
  load(ref)  # loads 'prototypes'
  names_to_keep <- prototypes$Accession
  
  # Require prototypes to be present in the input alignment (as documented)
  if (!all(names_to_keep %in% names(fastaData))) {
    stop(
      "To classify rhinovirus sequences, your input must include the prototype ",
      "sequences (same accessions). You have two options:\n",
      "  (1) Align in R with packaged references: ",
      "alignToRefs(seqData = fastaData, ...), then pass the result to assignTypes(); or\n",
      "  (2) Export prototypes with getPrototypeSeqs(), combine with your sequences, ",
      "align using your preferred tool, and re-import the curated alignment."
    )
  }
  
  # Compute pairwise distances
  distances <- pairwiseDistances(
    fastaData,
    model = model,
    deleteGapsGlobally = deleteGapsGlobally
  )
  
  # Keep ONLY prototype columns; remove prototype rows (so rows are queries)
  distances <- distances[, colnames(distances) %in% names_to_keep, drop = FALSE]
  distances <- distances[!rownames(distances) %in% names_to_keep, , drop = FALSE]
  
  # Assign one query (row) at a time
  assign_one <- function(row_vals) {
    # All NA row: nothing to decide
    if (all(is.na(row_vals))) {
      return(c("unassigned", NA_character_, NA_character_))
    }
    
    # Nearest prototype overall
    minCol <- which.min(row_vals)
    if (length(minCol) == 0L || is.na(row_vals[minCol])) {
      return(c("unassigned", NA_character_, NA_character_))
    }
    nearest_ref  <- colnames(distances)[minCol]
    nearest_dist <- row_vals[minCol]
    
    # Within threshold -> assign; otherwise keep unassigned but report distance
    if (!is.na(nearest_dist) && nearest_dist < threshold) {
      assignedType <- sub(".*_", "", gsub("RV", "", nearest_ref))
      return(c(assignedType, as.character(nearest_dist), nearest_ref))
    } else {
      return(c("unassigned", as.character(nearest_dist), nearest_ref))
    }
  }
  
  # Apply over all queries
  result <- t(apply(distances, 1, assign_one))
  
  # Assemble output
  outputDf <- data.frame(
    query        = rownames(distances),
    assignedType = result[, 1],
    distance     = as.numeric(result[, 2]),
    reference    = result[, 3],
    stringsAsFactors = FALSE
  )
  
  return(outputDf)
}
