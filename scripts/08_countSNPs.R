

source("scripts/04_genetic_distances.R")

countSNPs <- function(inputSequencesPath){
  pathToRef = inputSequencesPath
  queryFastaData = readFasta(inputSequencesPath)
  # run countSNP function
  snps <- countSNPsHelper(pathToRef, queryFastaData)
  # output
  return(snps)
}

# Example usage
countSNPs(inputSequencesPath = "./data/RVBPrototypeAligned.fasta")



