library(ape)
# read input
fastaData <- readFasta("./data/test.fasta", desiredLength = 420)

fastaData_ape <- read.dna('./data/test.fasta', format = 'fasta')


fastaData <- readFasta("./data/input_aln.fasta", desiredLength = 420)

fastaData_ape <- read.dna('./data/input_aln.fasta', format = 'fasta')

# assignTypes

#assigned <- assignTypes(fastaData)




# distances
dist_matrix1 <- pairwiseDistances(fastaData, model = "JC")
dist_matrix2_ape <- dist.dna(fastaData_ape, model = 'JC69')


heatmap(as.matrix(dist_matrix1), Rowv = NA, Colv = NA, scale = "none", 
        main = "Distance Matrix 1")
heatmap(as.matrix(dist_matrix2_ape), Rowv = NA, Colv = NA, scale = "none", 
        main = "Distance Matrix 2")



# Perform MDS analysis
mds1 <- cmdscale(as.dist(as.matrix(dist_matrix1)))
mds2 <- cmdscale(as.dist(as.matrix(dist_matrix2_ape)))


write.csv(as.matrix(dist_matrix2_ape), "~/Desktop/ape_output.csv")
write.csv(as.matrix(dist_matrix1), "~/Desktop/rhinotypeR_output.csv")


# Plot MDS
plot(mds1, type = 'n', main = "MDS Comparison")
points(mds1, col = 'red', pch = 19)
points(mds2, col = 'blue', pch = 18)
legend("topright", legend = c("Dataset 1", "Dataset 2"), 
       col = c("red", "blue"), pch = c(19, 18))





cor(as.vector(dist_matrix1), as.vector(dist_matrix2_ape))



# Ape
# Load the ape package
library(ape)

# Reading DNA sequences (assuming they are aligned)
# Replace 'file.fasta' with the path to your actual FASTA file
dna_sequences <- read.dna('./data/input_aln.fasta', format = 'fasta')

# Calculate pairwise distances
# 'model' can be one of several substitution models: "raw", "JC69", "K80", "F81", "K81", 
# "F84", "BH87", "T92", "TN93", "GG95", "logdet", "paralin", "indel"
pairwise_distances <- dist.dna(dna_sequences, model = 'JC69')

# View the pairwise distance matrix
pairwise_distances







# SeqinR
library(seqinr)

# Read in the alignment (mafft output)
dat <- read.alignment('data/input_aln.fasta',format = 'fasta')

# Calculate the pairwise distances
dist.alignment(dat,matrix='identity')