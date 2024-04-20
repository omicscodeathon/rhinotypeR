




source("scripts/04_genetic_distances.R")

pairwiseDistances <- function(fastaData, model = "p-distance") {

  # Determine which model to use based on user input
  if (model == "p-distance") {
    result <- calcPDistance(fastaData)
  } else if (model == "JC") {
    result <- calcJukesCantorDistance(fastaData)
  } else if (model == "Kimura2p") {
    result <- calcKimura2pDistance(fastaData)
  } else if (model == "Tamura3p") {
    result <- calcTamura3pDistance(fastaData)
  } else {
    stop("Unknown model specified. Choose from 'p-distance', 'JC', 'Kimura2p', or 'Tamura3p' ")
  }
  
  # Return the result of the chosen model
  return(result)
}


# Example usage
fastaD <- readFasta("./data/RVBPrototypeAligned.fasta")

pairwiseDistances(fastaD, "p-distance")
pairwiseDistances(fastaD, "JC")
pairwiseDistances(fastaD, "Kimura2p")
pairwiseDistances(fastaD, "Tamura3p")




