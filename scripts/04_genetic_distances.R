

# Base R functions

source("./scripts/02_readFasta.R")

# Function 4
# Count SNPs
countSNPsHelper <- function(fastaData) {
  
  # Read refs & query sequences to extract seq and header data 
  refs <- queries <- fastaData$sequences
  refHeaders <- queryHeaders <- fastaData$headers
  
  # Prepare the matrix to store SNP counts, using query headers for row names
  snpMatrix <- matrix(nrow = length(queries), ncol = length(refs), dimnames = list(queryHeaders, refHeaders))
  
  for (q in seq_along(queries)) {
    query <- queries[[q]]
    queryHeader <- queryHeaders[q]
    
    # Convert the query sequence to character vector
    queryChars <- strsplit(query, split = "")[[1]]
    
    for (i in seq_along(refs)) {
      # Convert the current reference sequence to character vector
      refChars <- strsplit(refs[[i]], split = "")[[1]]
      
      # Ensure both sequences are of the same length for comparison
      if (length(refChars) == length(queryChars)) {
        # Count the differences (SNPs)
        snps <- sum(refChars != queryChars)
      } else {
        snps <- NA  # Mark as NA if sequences are of different lengths
      }
      
      # Store SNP counts
      snpMatrix[q, i] <- snps
    }
  }
  
  # Output
  return(snpMatrix)
}


# Example usage
# Read in the query data first
queryFastaData <- readFasta("./data/tmp.fasta")
countSNPsHelper(queryFastaData)



# Function 5
# p-distance
calcPDistance <- function(fastaData) {
  # Count the SNPs between each query and each reference sequence
  snpCounts <- countSNPsHelper(fastaData)
  
  # Extract query headers
  queryHeaders <- refHeaders <- fastaData$headers
  
  # Directly calculate reference sequence lengths
  refLengths <- sapply(fastaData$sequences, function(seq) length(unlist(strsplit(seq, ""))))
  
  # Prepare a matrix for p-distances with appropriate dimensions and names
  pDistancesMatrix <- matrix(nrow = length(queryHeaders), ncol = length(refLengths), dimnames = list(queryHeaders, refHeaders))
  
  # Calculate p-distance for each query-reference pair
  for (q in 1:nrow(snpCounts)) {
    for (i in seq_along(refLengths)) {
      if (!is.na(snpCounts[q, i])) {
        pDistancesMatrix[q, i] <- snpCounts[q, i] / refLengths[i]
      } else {
        pDistancesMatrix[q, i] <- NA  # Assign NA if SNP count was NA
      }
    }
  }
  
  return(pDistancesMatrix)
}


# Example usage
# Read in the query data first
queryFastaData <- readFasta("./data/tmp.fasta")
calcPDistance(queryFastaData)


# Jukes Cantor

# Function to calculate Jukes-Cantor genetic distance
calcJukesCantorDistance <- function(fastaData) {
  # Calculate p-distance for multiple queries
  p_dist <- calcPDistance(fastaData)  # This now handles multiple queries
  
  # Initialize a matrix to store Jukes-Cantor distances
  jc_dist <- matrix(nrow = nrow(p_dist), ncol = ncol(p_dist), dimnames = dimnames(p_dist))
  
  # Apply the Jukes-Cantor formula to each element in the p-distance matrix
  jc_dist <- -3/4 * log(1 - 4/3 * p_dist)
  
  # Handling cases where p_dist >= 0.75, setting JC distance to Inf 
  # JC assumes that all nt substitutions are equally probable and independent, which might not always hold true in real data. 
      #The handling of cases where p_dist >= 0.75 with Inf is to indicate that the Jukes-Cantor model might not be valid for these high levels of divergence due to the assumption of the model being violated. 
  jc_dist[p_dist >= 0.75] <- Inf  
  
  # Return the Jukes-Cantor genetic distance matrix
  return(jc_dist)
}

# Example usage
# Read in the query data first
queryFastaData <- readFasta("./data/tmp.fasta")
calcJukesCantorDistance(queryFastaData)


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
calcKimura2pDistance <- function(fastaData) {
  
  # Load reference & query seq data and headers
  refs <- queries <- fastaData$sequences
  refHeaders <- queryHeaders <- fastaData$headers
  
  # Initialize a matrix to store the distances, with dimensions based on the number of queries and refs
  k2pMatrix <- matrix(nrow = length(queries), ncol = length(refs), 
                      dimnames = list(queryHeaders, refHeaders))
  
  # Iterate over each query sequence
  for (q in 1:length(queries)) {
    queryChars <- strsplit(queries[[q]], split = "")[[1]]
    
    # Iterate over each reference sequence
    for (r in 1:length(refs)) {
      refChars <- strsplit(refs[[r]], split = "")[[1]]
      
      # Ensure lengths match to avoid errors in distance calculation
      if (length(queryChars) == length(refChars)) {
        # Calculate proportions of transitions and transversions
        tt <- calcTransitionTransversions(refChars, queryChars)
        P <- tt$P
        Q <- tt$Q
        
        # Calculate Kimura 2-parameter genetic distance
        # Handling potential division by zero or negative values under the log function
        if ((1 - 2*P - Q) > 0 && (1 - 2*Q) > 0) {
          K2P_distance <- -0.5 * log((1 - 2*P - Q) * sqrt(1 - 2*Q))
        } else {
          K2P_distance <- NA  # Assign NA if the calculation is not possible
        }
        
        # Store the distance
        k2pMatrix[q, r] <- K2P_distance
      } else {
        k2pMatrix[q, r] <- NA  # Assign NA if sequence lengths do not match
      }
    }
  }
  
  # Return the matrix of Kimura 2-parameter genetic distances
  return(k2pMatrix)
}


# Example usage
# Read in the query data first
queryFastaData <- readFasta("./data/tmp.fasta")
calcKimura2pDistance(queryFastaData)


# calculate Tamura 3-parameter genetic distance
calcTamura3pDistance <- function(fastaData) {
  # Load reference & query seq data and headers
  refs <- queries <- fastaData$sequences
  refHeaders <- queryHeaders <- fastaData$headers
  
  # Initialize a matrix to store the distances, with dimensions based on the number of queries and refs
  tamuraMatrix <- matrix(nrow = length(queries), ncol = length(refs), 
                         dimnames = list(queryHeaders, refHeaders))
  
  # Iterate over each query sequence
  for (q in 1:length(queries)) {
    queryChars <- strsplit(queries[[q]], split = "")[[1]]
    
    # Iterate over each reference sequence
    for (r in 1:length(refs)) {
      refChars <- strsplit(refs[[r]], split = "")[[1]]
      
      # Ensure sequence lengths match for comparison
      if (length(queryChars) == length(refChars)) {
        # Calculate base frequencies and transitions/transversions
        gcContentRef <- sum(refChars %in% c("G", "C")) / length(refChars)
        gcContentQuery <- sum(queryChars %in% c("G", "C")) / length(queryChars)
        tt_results <- calcTransitionTransversions(refChars, queryChars)
        P <- tt_results$P
        Q <- tt_results$Q
        theta1 <- gcContentRef
        theta2 <- gcContentQuery
        C <- theta1 + theta2 - 2 * theta1 * theta2
        
        # Calculate Tamura 3-parameter distance with error handling
        if ((1 - P/C - Q) > 0 && (1 - 2*Q) > 0) {
          distance <- -C * log(1 - P/C - Q) - 0.5 * (1 - C) * log(1 - 2*Q)
        } else {
          distance <- NA  # Assign NA if distance cannot be calculated
        }
        
        tamuraMatrix[q, r] <- distance
      } else {
        tamuraMatrix[q, r] <- NA  # Assign NA if sequence lengths do not match
      }
    }
  }
  
  return(tamuraMatrix)
}

# Example usage
# Read in the query data first
queryFastaData <- readFasta("./data/tmp.fasta")
calcTamura3pDistance(queryFastaData)


# Tamura-Nei model  (complex!!!)

# General Time Reversible (GTR) Model (complex!!!)

# Maximum Composite Likelihood method? (possibly too complex for regular scripting!!!)





