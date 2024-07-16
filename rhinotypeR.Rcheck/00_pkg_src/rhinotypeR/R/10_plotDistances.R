source("R/05_pairwiseDistances.R")

plotDistances <- function(distancesMatrix){
  # Convert the data to a matrix 
  distancesMatrix <- as.matrix(distancesMatrix)
  
  # Heatmap
  heatmap(distancesMatrix, Rowv = NA, Colv = NA, col = heat.colors(256), 
          scale = "none", margins = c(7,10))
}

