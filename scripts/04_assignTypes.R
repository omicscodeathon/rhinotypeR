

source("./scripts/04_genetic_distances.R")

assignTypes <- function(pathToRef, pathToQuery, model = "p-distance", threshold = 0.105) {
  # Compute distances using the specified model
  distances <- allPrototypeDistances(pathToRef, pathToQuery, model)
  
  # Initialize vectors to store output data
  queryVec <- character(0)
  assignedTypeVec <- character(0)
  distanceVec <- numeric(0)
  
  # Iterate over each row (query) in the distances matrix
  for (i in 1:nrow(distances)) {
    # Extract row name (query header)
    queryHeader <- rownames(distances)[i]
    
    # Find columns (reference sequences) where distance is less than threshold
    validCols <- which(distances[i, ] < threshold)
    
    # If any valid columns found, add them to the vectors
    if (length(validCols) > 0) {
      for (col in validCols) {
        queryVec <- c(queryVec, queryHeader)
        assignedTypeVec <- c(assignedTypeVec, colnames(distances)[col])
        distanceVec <- c(distanceVec, distances[i, col])
      }
    }
  }
  
  # Create a data frame from the vectors
  outputDf <- data.frame(query = queryVec, assigned_type = assignedTypeVec, distance = distanceVec, stringsAsFactors = FALSE)
  
  # Write the data frame to a CSV file
  return(outputDf)
}

# Example usage
assignTypes("./data/RVBPrototypeAligned.fasta", "./data/tmp_query.fasta", "p-distance", 0.105)
