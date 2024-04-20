


source("./scripts/05_pairwiseDistances.R")


plotDistances <- function(distance_matrix){
  # Convert the data to a matrix 
  distance_matrix <- as.matrix(distance_matrix)
  
  # Heatmap
  heatmap(distance_matrix, Rowv = NA, Colv = NA, col = heat.colors(256), 
          scale = "none", margins = c(7,10))
  }



# Example usage
fastaD <- readFasta("./data/RVBPrototypeAligned.fasta")

distancesToPrototypes <- pairwiseDistances(fastaD, "p-distance")

plotDistances(distancesToPrototypes)
