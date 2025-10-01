#' Example VP4/2 alignment (DNAStringSet)
#'
#' A VP4/2 alignment used in examples and tests. Useful for demonstrating
#' `alignToRefs()`, `pairwiseDistances()`, and `assignTypes()`.
#'
#' @format A [Biostrings::DNAStringSet-class] where all sequences have equal
#'   width (aligned). Names are sequence IDs.
#'
#' @usage data(rhinovirusVP4)
#'
#' @details
#' This dataset is provided for illustration of the rhinotypeR workflow.
#' For related cohort data and background context on respiratory viruses in a
#' coastal Kenya school cohort, see the sources and references below.
#'
#' @examples
#' data(rhinovirusVP4)
#' rhinovirusVP4
#' # d <- pairwiseDistances(rhinovirusVP4, model = "IUPAC")
#'
#' @source
#' Luka, M. M., et al. (2020). Molecular Epidemiology of Human Rhinovirus 
#' From 1-Year Surveillance Within a School Setting in Rural Coastal Kenya.
#'
#' @references
#' Luka, M. M., et al. (2020). Molecular Epidemiology of Human Rhinovirus 
#' From 1-Year Surveillance Within a School Setting in Rural Coastal Kenya.
#' *Open Forum Infectious Diseases*, 7(10), ofaa385.
#' \doi{10.1093/ofid/ofaa385} \cr
#'
#' Adema, I. et al. (2020).Surveillance of respiratory viruses among children attending 
#' a primary school in rural coastal Kenya (Wellcome Open Research, 5:63, v2).
#' \doi{10.12688/wellcomeopenres.15703.2}\cr
#'
#' @keywords datasets
#' @name rhinovirusVP4
#' @docType data
"rhinovirusVP4"
