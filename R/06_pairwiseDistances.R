source("R/04_genetic_distances.R")

pairwiseDistances <- function(inputSequencesPath, model = "p-distance") {
  pathToRef = inputSequencesPath
  queryFastaData = readFasta(inputSequencesPath)
  
  # Determine which model to use based on user input
  if (model == "p-distance") {
    result <- calcPDistance(pathToRef, queryFastaData)
  } else if (model == "JC") {
    result <- calcJukesCantorDistance(pathToRef, queryFastaData)
  } else if (model == "Kimura2p") {
    result <- calcKimura2pDistance(pathToRef, queryFastaData)
  } else if (model == "Tamura3p") {
    result <- calcTamura3pDistance(pathToRef, queryFastaData)
  } else {
    stop("Unknown model specified. Choose from 'p-distance', 'JC', 'Kimura2p', or 'Tamura3p' ")
  }
  
  # Return the result of the chosen model
  return(result)
}


# Example usage
#pairwiseDistances("./data/RVBPrototypeAligned.fasta", "p-distance")
#pairwiseDistances("./data/RVBPrototypeAligned.fasta", "JC")
#pairwiseDistances("./data/RVBPrototypeAligned.fasta", "Kimura2p")
#pairwiseDistances("./data/RVBPrototypeAligned.fasta", "Tamura3p")




