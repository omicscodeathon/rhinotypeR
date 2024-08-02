## ----setup, include=FALSE-----------------------------------------------------
options(repos = c(CRAN = "https://cran.rstudio.com/"))


## ----eval=FALSE---------------------------------------------------------------
#  if (!require("BiocManager", quietly = TRUE))
#      install.packages("BiocManager")
#  
#  BiocManager::install("rhinotypeR")

## ----message=FALSE, warning=FALSE, results='hide'-----------------------------
install.packages("devtools", quietly = TRUE)

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

## ----echo=FALSE---------------------------------------------------------------
distances[1:5, 1:7]

## -----------------------------------------------------------------------------
genotypes <- assignTypes(sequences, model = "p-distance", gapDeletion = TRUE, threshold = 0.105)

head(genotypes)


## -----------------------------------------------------------------------------
plotFrequency(genotypes)

## -----------------------------------------------------------------------------
plotDistances(distances)

## -----------------------------------------------------------------------------
# sub-sample 
sampled_distances <- distances[1:30,1:30]

plotTree(sampled_distances)

