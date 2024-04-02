

source("scripts/04_genetic_distances.R")

countSNPs <- function(inputSequencesPath){
  pathToRef = inputSequencesPath
  pathToQuery = inputSequencesPath
  # run countSNP function
  snps <- countSNPsHelper(pathToRef, pathToQuery)
  # output
  return(snps)
}


countSNPs(inputSequencesPath = "./data/RVBPrototypeAligned.fasta")



