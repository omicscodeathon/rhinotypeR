
pairwiseDistances <- function(fastaData, model = "p-distance", gapDeletion = TRUE) {
  
  # preprocess fasta data
  fastaData <- preProcessFastaStringSet(fastaData)

  # Determine which model to use based on user input
  if (model == "p-distance") {
    result <- calcPDistance(fastaData, gapDeletion = gapDeletion)
  } else if (model == "JC") {
    result <- calcJukesCantorDistance(fastaData, gapDeletion = gapDeletion)
  } else if (model == "Kimura2p") {
    result <- calcKimura2pDistance(fastaData, gapDeletion = gapDeletion)
  } else if (model == "Tamura3p") {
    result <- calcTamura3pDistance(fastaData, gapDeletion = gapDeletion)
  } else {
    stop("Unknown model specified. Choose from 'p-distance', 'JC', 'Kimura2p', or 'Tamura3p' ")
  }
  
  # Return the result of the chosen model
  return(result)
}
