pkgname <- "rhinotypeR"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('rhinotypeR')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("SNPeek")
### * SNPeek

flush(stderr()); flush(stdout())

### Name: SNPeek
### Title: Visualize single nucleotide polymorphisms
### Aliases: SNPeek
### Keywords: SNP visualization

### ** Examples

# Load the dataset
test <- system.file("extdata", "test.fasta", package = "rhinotypeR")

fastaData <- readFasta(fastaFile = test)
SNPeek(fastaData, showLegend = FALSE)



cleanEx()
nameEx("assignTypes")
### * assignTypes

flush(stderr()); flush(stdout())

### Name: assignTypes
### Title: Assigns genotypes to input sequences
### Aliases: assignTypes
### Keywords: genotype sequence analysis

### ** Examples

# Load the dataset
test <- system.file("extdata", "input_aln.fasta", package = "rhinotypeR")

# Run command
fastaD <- readFasta(test)
assignTypes(fastaD, model = "p-distance", gapDeletion = TRUE, threshold = 0.105)



cleanEx()
nameEx("countSNPs")
### * countSNPs

flush(stderr()); flush(stdout())

### Name: countSNPs
### Title: Counts single nucleotide polymorphisms
### Aliases: countSNPs
### Keywords: genotype sequence analysis

### ** Examples

# Load the dataset
test <- system.file("extdata", "test.fasta", package = "rhinotypeR")

# Run the function
fastaData <- readFasta(test)
countSNPs(fastaData)



cleanEx()
nameEx("getPrototypeSeqs")
### * getPrototypeSeqs

flush(stderr()); flush(stdout())

### Name: getPrototypeSeqs
### Title: Download rhinovirus prototype strains
### Aliases: getPrototypeSeqs
### Keywords: genotype sequence analysis

### ** Examples

# Run the function
# getPrototypeSeqs(destinationFolder = "~/Desktop")



cleanEx()
nameEx("overallMeanDistance")
### * overallMeanDistance

flush(stderr()); flush(stdout())

### Name: overallMeanDistance
### Title: Estimates the overall mean distance
### Aliases: overallMeanDistance
### Keywords: genotype sequence analysis

### ** Examples

# Load the dataset
test <- system.file("extdata", "input_aln.fasta", package = "rhinotypeR")

# Usage
fastaData <- readFasta(test)
overallMeanDistance(fastaData, model="p-distance")



cleanEx()
nameEx("pairwiseDistances")
### * pairwiseDistances

flush(stderr()); flush(stdout())

### Name: pairwiseDistances
### Title: Estimates pairwise distances
### Aliases: pairwiseDistances
### Keywords: genotype sequence analysis

### ** Examples

# Load the dataset
test <- system.file("extdata", "input_aln.fasta", package = "rhinotypeR")

# Example usage
fastaD <- readFasta(test)
pairwiseDistances(fastaD, "p-distance")



cleanEx()
nameEx("plotAA")
### * plotAA

flush(stderr()); flush(stdout())

### Name: plotAA
### Title: Visualize amino acid substitutions
### Aliases: plotAA
### Keywords: genotype sequence analysis

### ** Examples

# Load the dataset
test <- system.file("extdata", "test.translated.fasta", package = "rhinotypeR")

# Usage
plotAA(test)



cleanEx()
nameEx("plotDistances")
### * plotDistances

flush(stderr()); flush(stdout())

### Name: plotDistances
### Title: Visualizes pairwise genetic distances
### Aliases: plotDistances
### Keywords: genotype sequence analysis

### ** Examples

# Load the dataset
test <- system.file("extdata", "input_aln.fasta", package = "rhinotypeR")

# Example usage
fastaD <- readFasta(test)
distancesMatrix <- pairwiseDistances(fastaD, "p-distance")
plotDistances(distancesMatrix)



cleanEx()
nameEx("plotFrequency")
### * plotFrequency

flush(stderr()); flush(stdout())

### Name: plotFrequency
### Title: Plots the frequency of assigned genotypes
### Aliases: plotFrequency
### Keywords: genotype visualization

### ** Examples

# Load the dataset
test <- system.file("extdata", "input_aln.fasta", package = "rhinotypeR")

# Run 
queryFastaData <- readFasta(test)
df <- assignTypes(queryFastaData, "p-distance")

plotFrequency(df)



cleanEx()
nameEx("plotTree")
### * plotTree

flush(stderr()); flush(stdout())

### Name: plotTree
### Title: Plots a simple phylogenetic tree
### Aliases: plotTree
### Keywords: genotype phylogenetics

### ** Examples

# Load the dataset
test <- system.file("extdata", "input_aln.fasta", package = "rhinotypeR")

# Example usage
fastaD <- readFasta(test)
pdistances <- pairwiseDistances(fastaD, "p-distance")
plotTree(pdistances)



cleanEx()
nameEx("readFasta")
### * readFasta

flush(stderr()); flush(stdout())

### Name: readFasta
### Title: Reads sequence alignment/fasta files into R for processing
### Aliases: readFasta
### Keywords: sequence analysis data input

### ** Examples

# Load the dataset
test <- system.file("extdata", "test.fasta", package = "rhinotypeR")

# Run the command
readFasta(test)



### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
