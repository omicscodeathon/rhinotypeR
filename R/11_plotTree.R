plotTree <- function(distance_matrix, ...) {
  
  # Convert the data to a matrix 
  distance_matrix <- as.matrix(distance_matrix)
  
  # Convert the distance matrix into a "dist" object required by hclust
  distance_object <- as.dist(distance_matrix)
  
  # Perform hierarchical clustering using complete linkage
  hc <- hclust(distance_object, method = "complete")
  
  # Plot the dendrogram with additional arguments passed through `...`
  plot(hc, ...)
}
