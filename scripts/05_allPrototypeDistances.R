


source("scripts/03_genetic_distances.R")

allPrototypeDistances <- function(pathToRef, pathToQuery, model = "p-distance") {
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
    stop("Unknown model specified. Choose from 'p-distance', 'JC', 'Kimura2p', or 'Tamura3p'.")
  }
  
  # Return the result of the chosen model
  return(result)
}



# Example usage
allPrototypeDistances("./data/RVBPrototypeAligned.fasta", "./data/tmp_query.fasta", "p-distance")
allPrototypeDistances("./data/RVBPrototypeAligned.fasta", "./data/tmp_query.fasta", "JC")
allPrototypeDistances("./data/RVBPrototypeAligned.fasta", "./data/tmp_query.fasta", "Kimura2p")
allPrototypeDistances("./data/RVBPrototypeAligned.fasta", "./data/tmp_query.fasta", "Tamura3p")

