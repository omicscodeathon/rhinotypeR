

countSNPs <- function(fastaData, gapDeletion = TRUE){
  # preprocess fasta data
  fastaData <- preProcessFastaStringSet(fastaData)
  
  # run countSNP function
  snps <- countSNPsHelper(fastaData, gapDeletion = gapDeletion)
  # output
  return(snps)
}

