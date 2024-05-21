


source("./scripts/06_pairwiseDistances.R")


plotTree <- function(distance_matrix) {
  
  # Convert the data to a matrix 
  distance_matrix <- as.matrix(distance_matrix)
  
  # Convert the distance matrix into a "dist" object required by hclust
  distance_object <- as.dist(distance_matrix)
  
  # Perform hierarchical clustering using complete linkage
  hc <- hclust(distance_object, method = "complete")
  
  # Plot the dendrogram
  plot(hc, hang = -1, cex = 0.6, main = "A simple tree", xlab = "", #ann = par("ann"),
       ylab = "Genetic distance")
}


# Example usage
fastaD <- readFasta("./data/input_aln.fasta")

pdistances <- pairwiseDistances(fastaD, "p-distance")

plotTree(pdistances)

