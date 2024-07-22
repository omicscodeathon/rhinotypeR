## ----setup, include=FALSE-----------------------------------------------------
options(repos = c(CRAN = "https://cran.rstudio.com/"))


## ----eval=FALSE---------------------------------------------------------------
#  if (!require("BiocManager", quietly = TRUE))
#      install.packages("BiocManager")
#  
#  BiocManager::install("Biostrings")

## -----------------------------------------------------------------------------
install.packages("devtools")

devtools::install_github("omicscodeathon/rhinotypeR")

## -----------------------------------------------------------------------------
library(rhinotypeR)

## ----eval=FALSE---------------------------------------------------------------
#  getPrototypeSeqs("~/Desktop")

## -----------------------------------------------------------------------------
sequences <- readFasta("../inst/extdata/input_aln.fasta")

## -----------------------------------------------------------------------------
SNPeek(sequences)

## -----------------------------------------------------------------------------
distances <- pairwiseDistances(sequences, model = "p-distance", gapDeletion = TRUE)

## -----------------------------------------------------------------------------
genotypes <- assignTypes(sequences, model = "p-distance", gapDeletion = TRUE, threshold = 0.105)

## -----------------------------------------------------------------------------
plotFrequency(genotypes)
plotDistances(distances)

## -----------------------------------------------------------------------------
plotDistances(distances)

## -----------------------------------------------------------------------------
plotTree(distances)

