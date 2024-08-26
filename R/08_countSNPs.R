

countSNPs <- function(fastaData, gapDeletion = TRUE){
  # run countSNP function
  snps <- countSNPsHelper(fastaData, gapDeletion = gapDeletion)
  # output
  return(snps)
}

