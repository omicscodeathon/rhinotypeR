# Calculate pairwise distances

##Load libraries
library(ape)
library(seqinr)

##read in data
dat <- read.alignment('data/Ref_target.fasta',format = 'fasta')

# -------------------------------------------------------------------------
##Method 1 - dist.alignment
dist.alignment(dat,matrix='identity')

# -------------------------------------------------------------------------
##Method 2 - pairdist_seqinr
pairdist_seqinr <- function(data,matrix){
  library(seqinr)
  pd=dist.alignment(data,matrix)
  return(pd)
}
pairdist_seqinr(data=dat,matrix="identity")

# -------------------------------------------------------------------------
##Method 3 using ape
mat <- as.matrix(dat) #ape takes in a matrix, hence convert the alignment to a matrix
ape::dist.gene(mat,method='pairwise')

pairdist_ape <- function(data,method){
  library(ape)
  mat= as.matrix(data)
  pdape=ape::dist.gene(mat, method)
  return(pdape)
}

pairdist_ape(data=dat, method='pairwise')

# -------------------------------------------------------------------------
#Method 4 using dist.dna

# Extract DNA sequences from the read data
dna_data <- as.DNAbin(dat)

# Compute distances between DNA sequences
distances <- dist.dna(dna_data)

# Print the distance matrix
print(distances)


