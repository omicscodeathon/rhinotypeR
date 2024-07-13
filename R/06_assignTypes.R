source("R/04_genetic_distances.R")
source("R/05_pairwiseDistances.R")

assignTypes <- function(fastaData, model = "p-distance", gapDeletion = TRUE, threshold = 0.105) {
  
  ref <- system.file("extdata", "prototypes.csv", package = "rhinotypeR")
  prototypes <- read.csv(ref)
  names_to_keep <- prototypes$Accession
  
  # run pairwiseDistances to calculate distances
  distances <- pairwiseDistances(fastaData, model = model, gapDeletion=gapDeletion)
  
  # Filter columns based on the prototypes
  distances <- distances[, colnames(distances) %in% names_to_keep]
  
  # Filter rows to remove the same names as in the list
  distances <- distances[!rownames(distances) %in% names_to_keep, ]
  # Initialize vectors to store output data
  queryVec <- character(0)
  assignedTypeVec <- character(0)
  distanceVec <- numeric(0)
  refSeqVec <- character(0)
  
  # Iterate over each row (query) in the distances matrix
  for (i in seq_len(nrow(distances))) {
    queryHeader <- rownames(distances)[i]
    validCols <- which(distances[i, ] < threshold)
    
    if (length(validCols) == 0) {
      # If no valid columns found, mark as "unassigned"
      assignedTypeVec <- c(assignedTypeVec, "unassigned")
      distanceVec <- c(distanceVec, NA)
      # Find the column with the minimum distance
      minDistCol <- which.min(distances[i, ])
      refSeqVec <- c(refSeqVec, colnames(distances)[minDistCol])
    } else {
      # Choose the one with the minimum distance
      minDistanceCol <- which.min(distances[i, validCols])
      col <- validCols[minDistanceCol]
      assignedType <- colnames(distances)[col]
      assignedTypeVec <- c(assignedTypeVec, sub(".*_", "", gsub("RV", "", assignedType)))
      distanceVec <- c(distanceVec, distances[i, col])
      refSeqVec <- c(refSeqVec, assignedType)
    }
    queryVec <- c(queryVec, queryHeader)
  }
  
  outputDf <- data.frame(query = queryVec, assignedType = assignedTypeVec, 
                         distance = distanceVec, reference = refSeqVec, stringsAsFactors = FALSE)
  
  return(outputDf)
}

