


source("./scripts/06_pairwiseDistances.R")


distancesToPrototypes <- pairwiseDistances("./data/RVBPrototypeAligned.fasta", 
                                                "p-distance")


# Convert the data to a matrix 
distance_matrix <- as.matrix(distancesToPrototypes)

# Heatmap
heatmap(distance_matrix, Rowv = NA, Colv = NA, col = heat.colors(256), 
        scale = "none", margins = c(7,10))
