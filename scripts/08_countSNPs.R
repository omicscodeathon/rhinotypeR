

source("scripts/04_genetic_distances.R")

countSNPs <- function(fastaData, gapDeletion = TRUE){
  # run countSNP function
  snps <- countSNPsHelper(fastaData, gapDeletion = gapDeletion)
  # output
  return(snps)
}

# Example usage
source("./scripts/02_readFasta.R")
fastaData <- readFasta("./data/test.fasta")
countSNPs(fastaData)



