




source("scripts/04_genetic_distances.R")

queryPairwiseDistances <- function(inputSequencesPath, model = "p-distance") {
  pathToRef = inputSequencesPath
  pathToQuery = inputSequencesPath
  
  # Determine which model to use based on user input
  if (model == "p-distance") {
    result <- calcPDistance(pathToRef, pathToQuery)
  } else if (model == "JC") {
    result <- calcJukesCantorDistance(pathToRef, pathToQuery)
  } else if (model == "Kimura2p") {
    result <- calcKimura2pDistance(pathToRef, pathToQuery)
  } else if (model == "Tamura3p") {
    result <- calcTamura3pDistance(pathToRef, pathToQuery)
  } else {
    stop("Unknown model specified. Choose from 'p-distance', 'JC', 'Kimura2p', or 'Tamura3p' ")
  }
  
  # Return the result of the chosen model
  return(result)
}



# Example usage
queryPairwiseDistances("./data/RVBPrototypeAligned.fasta", "p-distance")
queryPairwiseDistances("./data/RVBPrototypeAligned.fasta", "JC")
queryPairwiseDistances("./data/RVBPrototypeAligned.fasta", "Kimura2p")
queryPairwiseDistances("./data/RVBPrototypeAligned.fasta", "Tamura3p")




