source("R/04_genetic_distances.R")

countSNPs <- function(fastaData){
  # run countSNP function
  snps <- countSNPsHelper(fastaData)
  # output
  return(snps)
}