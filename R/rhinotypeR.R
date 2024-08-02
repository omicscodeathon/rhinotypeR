#' rhinotypeR: A Package for the genotyping of rhinoviruses using the VP4/2 region
#'
#' The \pkg{rhinotypeR} package provides tools and functions that streamline the genotyping of rhinoviruses using the VP4/2 region.
#'
#' This package includes functions to read sequence data, calculate pairwise distances and visualization, specifically tailored for rhinovirus data. It also provides additional utilities for working with sequence data in FASTA format such as simple phylogenetic trees.
#'
#' @section Main functions:
#' \itemize{
#'   \item \code{\link{getPrototypeSeqs}}: Download prototype sequences into local machine.
#'   \item \code{\link{readFasta}}: Read sequences from a FASTA file.
#'   \item \code{\link{pairwiseDistances}}: Calculates pairwise distances using a specified evolutionary model.
#'   \item \code{\link{assignTypes}}: Assigns sequences to their respective rhinovirus genotypes.
#'   \item \code{\link{plotTree}}: Plots a simple phylogenetic tree.
#' }
#'
#' @section Getting started:
#' To get started with \pkg{rhinotypeR}, download the prototype sequences and combine these with your newly generated VP4/2 sequences and align using a suitable tool. Then import the curated alignment into R:
#' \preformatted{
#'
#' \# Download prototype sequences
#' getPrototypeSeqs("path/to/destination")
#'
#' \# Read sequences from a FASTA file
#' sequences <- readFasta("../inst/extdata/input_aln.fasta")
#'
#' \# Perform sequence analysis
#' distance_matrix <- pairwiseDistances(sequences)
#'
#' \# Assign to genotypes
#' assigned_types <- assignTypes(distance_matrix)
#'
#' \# Simple phylogenetic tree
#' plotTree(distance_matrix)
#' }
#'
#' For more detailed examples and usage, refer to the package vignettes:
#' \preformatted{
#' vignette(package = "rhinotypeR")
#' vignette("rhinotypeR", package = "rhinotypeR")
#' }
#'
#' @docType package
#' @name rhinotypeR
#' @aliases rhinotypeR-package
#' @seealso \url{https://github.com/omicscodeathon/rhinotypeR}
#' @keywords package
"_PACKAGE"
