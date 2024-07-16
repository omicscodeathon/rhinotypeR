## -----------------------------------------------------------------------------
# If needed
# install.packages("devtools")
# devtools::install_github("yourusername/rhinotypeR")

## -----------------------------------------------------------------------------
library(rhinotypeR)

## -----------------------------------------------------------------------------
# getPrototypeSeqs("path_to_destination_folder")
# For example
# getPrototypeSeqs("~/Desktop")

## -----------------------------------------------------------------------------
# sequences <- readFasta("path/to/your/fasta/file")

## -----------------------------------------------------------------------------
# SNPeek(sequences)

## -----------------------------------------------------------------------------
# distances <- pairwiseDistances(sequences)

## -----------------------------------------------------------------------------
# genotypes <- assignTypes(sequences)

## -----------------------------------------------------------------------------
# plotFrequency(genotypes)
# plotDistances(distances)

## -----------------------------------------------------------------------------
sessionInfo()

