#' Rhinovirus VP4/2 prototype references (DNAStringSet)
#'
#' A FASTA bundle of rhinovirus prototype sequences used as references for
#' genotyping (VP4/2 region). Provided as a convenience dataset for examples,
#' vignettes, and tests.
#'
#' @format A [Biostrings::DNAStringSet-class] containing near-aligned
#'   VP4/2 prototype sequences. Sequence names are accessions (optionally with
#'   type labels).
#'
#' @details
#' These sequences mirror the prototypes shipped in \file{inst/extdata/prototypes.fasta}.
#' You may use this object directly in workflows (e.g., to append to user sequences
#' before alignment) or export prototypes to disk with [getPrototypeSeqs()].
#'
#' @usage data(rhinovirusPrototypesVP4)
#'
#' @examples
#' data(rhinovirusPrototypesVP4)
#' rhinovirusPrototypesVP4
#'
#' @source McIntyre, C. L., Knowles, N. J., & Simmonds, P. (2013).
#' Proposals for the classification of human rhinovirus species A, B and C
#' into genotypically assigned types.
#'
#' @references
#' McIntyre, C. L., Knowles, N. J., & Simmonds, P. (2013).
#' Proposals for the classification of human rhinovirus species A, B and C
#' into genotypically assigned types. *Journal of General Virology*, 94(8),
#' 1791–1806. \doi{10.1099/vir.0.053686-0}
#'
#' @keywords datasets
#' @name rhinovirusPrototypesVP4
#' @docType data
"rhinovirusPrototypesVP4"

