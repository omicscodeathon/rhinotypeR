

source("./scripts/04_genetic_distances.R")

      #. Choose from 'p-distance', 'JC', 'Kimura2p', or 'Tamura3p'

assignTypes <- function(pathToRef, pathToQuery, model = "p-distance", threshold = 0.105) {
  # run allPrototypeDistances to calculate distances
  distances <- allPrototypeDistances(pathToRef, pathToQuery, model)
  
  # Initialize vectors to store output data
  queryVec <- character(0)
  assignedTypeVec <- character(0)
  distanceVec <- numeric(0)
  
  # Iterate over each row (query) in the distances matrix
  for (i in 1:nrow(distances)) {
    queryHeader <- rownames(distances)[i]
    
    validCols <- which(distances[i, ] < threshold)
    
    if (length(validCols) == 0) {
      # If no valid columns found, mark as "unassigned"
      queryVec <- c(queryVec, queryHeader)
      assignedTypeVec <- c(assignedTypeVec, "unassigned")
      distanceVec <- c(distanceVec, NA)
    } else if (length(validCols) > 1) {
      # If multiple valid columns, choose the one with the minimum distance
      minDistanceCol <- which.min(distances[i, validCols])
      col <- validCols[minDistanceCol]
      
      queryVec <- c(queryVec, queryHeader)
      assignedType <- colnames(distances)[col]
      assignedTypeCleaned <- sub(".*_", "", assignedType)
      assignedTypeCleaned <- gsub("RV", "", assignedTypeCleaned)
      
      assignedTypeVec <- c(assignedTypeVec, assignedTypeCleaned)
      distanceVec <- c(distanceVec, distances[i, col])
    } else {
      # For a single valid column
      col <- validCols
      
      queryVec <- c(queryVec, queryHeader)
      assignedType <- colnames(distances)[col]
      assignedTypeCleaned <- sub(".*_", "", assignedType)
      assignedTypeCleaned <- gsub("RV", "", assignedTypeCleaned)
      
      assignedTypeVec <- c(assignedTypeVec, assignedTypeCleaned)
      distanceVec <- c(distanceVec, distances[i, col])
    }
  }
  
  outputDf <- data.frame(query = queryVec, assigned_type = assignedTypeVec, distance = distanceVec, stringsAsFactors = FALSE)
  
  return(outputDf)
}

# Example usage
assignTypes("./data/RVBPrototypeAligned.fasta", "./data/tmp.fasta", "p-distance", 0.105)
assignTypes("./data/RVBPrototypeAligned.fasta", "./data/tmp_off_sequence.fasta", "Tamura3p", 0.105)
