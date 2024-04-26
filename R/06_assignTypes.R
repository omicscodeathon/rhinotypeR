 source("R/04_genetic_distances.R")
 source("R/05_pairwiseDistances.R")

#. Choose from 'p-distance', 'JC', 'Kimura2p', or 'Tamura3p'

assignTypes <- function(fastaData, model = "p-distance", threshold = 0.105) {
  
  # Load the system dataset
  ref <- system.file("extdata", "prototypes.csv", package = "rhinotypeR")
  
  prototypes <- read.csv(ref)
  
  names_to_keep <- prototypes$Accession
  
  # run pairwiseDistances to calculate distances
  distances <- pairwiseDistances(fastaData, model = model)
  
  # Filter columns based on the prototypes
  distances <- distances[, colnames(distances) %in% names_to_keep]
  
  # Filter rows to remove the same names as in the list
  distances <- distances[!rownames(distances) %in% names_to_keep, ]
  
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
