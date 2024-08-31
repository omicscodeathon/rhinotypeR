assignTypes <- function(fastaData, model = "p-distance", gapDeletion = TRUE, threshold = 0.105) {
  
  # Preprocess fasta data
  fastaData <- preProcessFastaStringSet(fastaData)
  
  # Load prototype sequences
  ref <- system.file("extdata", "prototypes.csv", package = "rhinotypeR")
  prototypes <- read.csv(ref)
  names_to_keep <- prototypes$Accession
  
  # Check if input fastaData contains the prototype sequences
  if (!all(names_to_keep %in% fastaData$headers)) {
    stop("To classify rhinovirus sequences, please ensure your input sequences contain the prototypes. 
         These can be downloaded using `getPrototypeSeqs`.")
  }
  
  # Run pairwiseDistances to calculate distances
  distances <- pairwiseDistances(fastaData, model = model, gapDeletion = gapDeletion)
  
  # Filter columns and rows based on the prototypes
  distances <- distances[, colnames(distances) %in% names_to_keep, drop = FALSE]
  distances <- distances[!rownames(distances) %in% names_to_keep, , drop = FALSE]
  
  # Function to assign a type based on distance for a single query
  assign_type <- function(row) {
    validCols <- which(row < threshold)
    
    if (length(validCols) == 0) {
      minDistCol <- which.min(row)
      return(c("unassigned", NA, colnames(distances)[minDistCol]))
    } else {
      minDistanceCol <- which.min(row[validCols])
      col <- validCols[minDistanceCol]
      assignedType <- colnames(distances)[col]
      return(c(sub(".*_", "", gsub("RV", "", assignedType)), row[col], assignedType))
    }
  }
  
  # Apply the assign_type function over each row
  result <- t(apply(distances, 1, assign_type))
  
  # Convert result to a data frame
  outputDf <- data.frame(query = rownames(distances), 
                         assignedType = result[, 1], 
                         distance = as.numeric(result[, 2]), 
                         reference = result[, 3], 
                         stringsAsFactors = FALSE)
  
  return(outputDf)
}
