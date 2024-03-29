source("R/genetic_distances.R")

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


# -------------------------------------------------------------------------
PlotFrequency <- function(genotypeassigned, xlab='The assigned RV genotype',ylab='Sample genotype frequency'){
  # This function takes assignTypes() output and the user-supplied x and y-axis labels.
  # It tabulates the frequency of unique entries of the assigned_type column
  # Visualize using a bar plot the frequency of each genotype present in user samples
  plot=genotypeassigned %>%
    count(assigned_type) %>%
    ggplot(aes(x=assigned_type,y=n,fill=assigned_type))+
    geom_col()+
    labs(x=xlab, y=ylab)+
    theme(axis.text.x = element_text(angle = 30, hjust = 0.5, vjust = 0.5))
  return(plot)
}


#Save data in rda format
#pathToRef <- "./data/RVBPrototypeAligned.fasta"
#usethis::use_data(pathToRef)
#pathToQuery <- "./data/tmp_query.fasta"
#usethis::use_data(pathToQuery)
#genotypeassigned <- assignTypes(pathToRef, pathToQuery, model = "p-distance", threshold = 0.105)
#usethis::use_data(genotypeassigned)


# Example usage
#assignTypes(pathToRef, pathToQuery, "p-distance", 0.105)
