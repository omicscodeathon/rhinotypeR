source("R/04_genetic_distances.R")

pairwiseDistances <- function(fastaData, model = "p-distance", pairwiseDeletion = FALSE) {
  
  # Typically, Kimura2p and Tamura3p dont consider gaps
  if (model == "Kimura2p" | model == "Tamura3p" ){
    pairwiseDeletion = FALSE
  }else{
    pairwiseDeletion = pairwiseDeletion
  }
  
  # Determine which model to use based on user input
  if (model == "p-distance") {
    result <- calcPDistance(fastaData, pairwiseDeletion)
  } else if (model == "JC") {
    result <- calcJukesCantorDistance(fastaData, pairwiseDeletion)
  } else if (model == "Kimura2p") {
    result <- calcKimura2pDistance(fastaData, pairwiseDeletion)
  } else if (model == "Tamura3p") {
    result <- calcTamura3pDistance(fastaData, pairwiseDeletion)
  } else {
    stop("Unknown model specified. Choose from 'p-distance', 'JC', 'Kimura2p', or 'Tamura3p' ")
  }
  
  # Return the result of the chosen model
  return(result)
}
