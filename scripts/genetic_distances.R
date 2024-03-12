

# Base R functions


# Function 1 
# Read fasta
readFasta <- function(fastaFile) {
  lines <- readLines(fastaFile)
  sequence <- paste(lines[-1], collapse = "")
  return(toupper(sequence)) # Convert sequence to uppercase
}

# Example usage
readFasta("./data/tmp_query.fasta")
readFasta("./data/tmp_ref_b99.fasta")

# Function 2
# Compare sequence lengths
compareLengths <- function(pathToRef, pathToQuery){
  ref <- readFasta(pathToRef)
  query <- readFasta(pathToQuery)
  
  # Convert the sequences to character vectors
  refLength <- length(unlist(strsplit(ref, "")))
  queryLength <- length(unlist(strsplit(query, "")))
  
  # Return both lengths for further use
  return(list(refLength = refLength, queryLength = queryLength))
}
# Example usage
compareLengths(pathToRef = "./data/tmp_ref_b99.fasta", pathToQuery="./data/tmp_query.fasta")

# Function 3
# Function to handle gaps or non-ATCG characters in sequences
handleGaps <- function(refSeq, querySeq, option = "delete") {
  # Convert sequences to character vectors
  refChars <- strsplit(toupper(refSeq), "")[[1]]
  queryChars <- strsplit(toupper(querySeq), "")[[1]]
  
  # Identify positions with non-ATCG characters
  valid_positions <- which((refChars %in% c("A", "T", "C", "G")) & (queryChars %in% c("A", "T", "C", "G")))
  
  if (option == "delete") {
    # Keep only positions with ATCG characters in both sequences
    refChars <- refChars[valid_positions]
    queryChars <- queryChars[valid_positions]
  } # If option is "ignore", we do nothing to the sequences
  
  # Reassemble sequences
  refSeq_cleaned <- paste(refChars, collapse = "")
  querySeq_cleaned <- paste(queryChars, collapse = "")
  
  return(list(refSeq = refSeq_cleaned, querySeq = querySeq_cleaned))
}



# Function 4
# Count SNPs
countSNPs <- function(pathToRef, pathToQuery){
  
  # Read fasta file
  ref <- readFasta(pathToRef)
  query <- readFasta(pathToQuery)
  
  # Convert the sequences to character vectors
  refChars <- strsplit(ref, split = "")[[1]]
  queryChars <- strsplit(query, split = "")[[1]]
  
  # Count the differences (SNPs)
  snps <- sum(refChars != queryChars)
  
  # output
  return(snps)
}
# Example usage
countSNPs(pathToRef = "./data/tmp_ref_b99.fasta", pathToQuery="./data/tmp_query.fasta")




# Function 5
# p-distance
calcPDistance <- function(pathToRef, pathToQuery) {
  
  # Count the differences (SNPs)
  snps <- countSNPs(pathToRef, pathToQuery)
  
  # lengths
  lengths <- compareLengths(pathToRef, pathToQuery)
  refLength <- lengths$refLength
  
  # Calculate p-distance
  p_dist <- snps/refLength
  
  # output
  return(p_dist)
}

# Example usage
calcPDistance(pathToRef = "./data/tmp_ref_b99.fasta", pathToQuery = "./data/tmp_query.fasta")


# Jukes Cantor

# Function to calculate Jukes-Cantor genetic distance
calcJukesCantorDistance <- function(pathToRef, pathToQuery) {
  # Calculate p-distance
  p_dist <- calcPDistance(pathToRef, pathToQuery)
  
  # Calculate Jukes-Cantor genetic distance
  jc_dist <- -3/4 * log(1 - 4/3 * p_dist)
  
  # Return the Jukes-Cantor genetic distance
  return(jc_dist)
}

# Example usage
calcJukesCantorDistance(pathToRef = "./data/tmp_ref_b99.fasta", pathToQuery = "./data/tmp_query.fasta")


# Kimura 2 parameter (transitions vs transversions)

## Function to calculate proportions of transitions and transversions
calcTransitionTransversions <- function(refChars, queryChars) {
  transitions <- sum((refChars == 'A' & queryChars == 'G') |
                       (refChars == 'G' & queryChars == 'A') |
                       (refChars == 'C' & queryChars == 'T') |
                       (refChars == 'T' & queryChars == 'C'))
  
  transversions <- sum((refChars == 'A' & queryChars == 'C') |
                         (refChars == 'A' & queryChars == 'T') |
                         (refChars == 'G' & queryChars == 'C') |
                         (refChars == 'G' & queryChars == 'T') |
                         (refChars == 'C' & queryChars == 'A') |
                         (refChars == 'C' & queryChars == 'G') |
                         (refChars == 'T' & queryChars == 'A') |
                         (refChars == 'T' & queryChars == 'G'))
  
  total_sites <- length(refChars)
  P <- transitions / total_sites
  Q <- transversions / total_sites
  
  return(list(P = P, Q = Q))
}

## Function to calculate Kimura 2-parameter genetic distance
calcKimura2pDistance <- function(pathToRef, pathToQuery) {
  
  # Read fasta files and ensure sequences are in uppercase
  ref <- readFasta(pathToRef)
  query <- readFasta(pathToQuery)
  
  # Convert the sequences to character vectors
  refChars <- strsplit(ref, split = "")[[1]]
  queryChars <- strsplit(query, split = "")[[1]]
  
  # Calculate proportions of transitions and transversions
  tt <- calcTransitionTransversions(refChars, queryChars)
  P <- tt$P
  Q <- tt$Q
  
  # Calculate Kimura 2-parameter genetic distance
  K2P_distance <- -0.5 * log((1 - 2*P - Q) * sqrt(1 - 2*Q))
  
  # Return the Kimura 2-parameter genetic distance
  return(K2P_distance)
}

# Example usage
calcKimura2pDistance(pathToRef = "./data/tmp_ref_b99.fasta", pathToQuery = "./data/tmp_query.fasta")


# calculate Tamura 3-parameter genetic distance
calcTamura3pDistance <- function(pathToRef, pathToQuery) {
  # Assuming readFasta is defined elsewhere
  ref <- toupper(readFasta(pathToRef))
  query <- toupper(readFasta(pathToQuery))
  
  # Convert sequences to character vectors
  refChars <- strsplit(ref, "")[[1]]
  queryChars <- strsplit(query, "")[[1]]
  
  # Calculate base frequencies for both sequences
  gcContentRef <- sum(refChars %in% c("G", "C")) / length(refChars)
  gcContentQuery <- sum(queryChars %in% c("G", "C")) / length(queryChars)
  
  # Use calcTransitionTransversions to calculate transitions and transversions
  tt_results <- calcTransitionTransversions(refChars, queryChars)
  P <- tt_results$P
  Q <- tt_results$Q
  
  theta1 <- gcContentRef
  theta2 <- gcContentQuery
  C <- theta1 + theta2 - 2 * theta1 * theta2
  
  # Tamura distance calculation
  if ((1 - P/C - Q) > 0 && (1 - 2*Q) > 0) {
    distance <- -C * log(1 - P/C - Q) - 0.5 * (1 - C) * log(1 - 2*Q)
  } else {
    distance <- NA  # Distance cannot be calculated under these conditions
  }
  
  return(distance)
}


# Example usage
calcTamura3pDistance(pathToRef = "./data/tmp_ref_b99.fasta", pathToQuery = "./data/tmp_query.fasta")


# Tamura-Nei model  (complex!!!)

# General Time Reversible (GTR) Model (complex!!!)

# Maximum Composite Likelihood method? (possibly too complex for regular scripting!!!)




assignTypes <- function(pathToRef, pathToQuery, model="p-distance"){
  
  
  
}



getPrototypeSeqs <- function(species="A", dest_path){
  
  
  
}


overallMeanDistance <- function(queryMultipleSequenceAlignment){
  
  
  
}


