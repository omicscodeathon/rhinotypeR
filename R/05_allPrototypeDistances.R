source("R/04_genetic_distances.R")

allPrototypeDistances <- function(pathToRef, queryFastaData, model = "p-distance") {
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
#allPrototypeDistances("./data/RVBPrototypeAligned.fasta", queryFastaData, "p-distance")
#allPrototypeDistances("./data/RVBPrototypeAligned.fasta", queryFastaData, "JC")
#allPrototypeDistances("./data/RVBPrototypeAligned.fasta", queryFastaData, "Kimura2p")
#allPrototypeDistances("./data/RVBPrototypeAligned.fasta", queryFastaData, "Tamura3p")

