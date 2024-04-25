

# read input
fastaData <- readFasta("./data/input_aln.fasta")

# assignTypes

assigned <- assignTypes(fastaData)







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



