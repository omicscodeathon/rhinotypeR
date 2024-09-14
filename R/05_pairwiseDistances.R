pairwiseDistances <- function(fastaData, model = "p-distance", gapDeletion = TRUE) {
  
  if (inherits(fastaData, "DNAMultipleAlignment") || inherits(fastaData, "DNAStringSet")) {
    # If it is either a DNAMultipleAlignment or DNAStringSet, preprocess the fasta data
    fastaData <- preProcessFastaStringSet(fastaData)
  }
  
  # Map model names to their corresponding functions
  functionMap <- list(
    "p-distance" = calcPDistance,
    "JC" = calcJukesCantorDistance,
    "Kimura2p" = calcKimura2pDistance,
    "Tamura3p" = calcTamura3pDistance
  )
  
  # Apply the appropriate function based on the model
  result <- applyModelFunction(fastaData, model, gapDeletion, functionMap)
  
  return(result)
}
