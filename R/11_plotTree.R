source("R/06_pairwiseDistances.R")

plotTree <- function(inputSequencesPath, model = "p-distance") {
  
  distancesToPrototypes <- pairwiseDistances(inputSequencesPath, "p-distance")
  
  # Convert the data to a matrix 
  distance_matrix <- as.matrix(distancesToPrototypes)
  
  # Convert the distance matrix into a "dist" object which is required by hclust
  distance_object <- as.dist(distance_matrix)
  
  # Perform hierarchical clustering using complete linkage
  hc <- hclust(distance_object, method = "complete")
  
  # Plot the dendrogram
  plot(hc, hang = -1, cex = 0.6, main = "A simple tree", xlab = "", #ann = par("ann"),
       ylab = "Genetic distance")

}




