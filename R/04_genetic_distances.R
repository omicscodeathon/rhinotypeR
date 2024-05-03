source("R/02_readFasta.R")


# Optionally delete sites with at least one missing data for all sequences 
deleteMissingDataSites <- function(seqs) {
  # Convert list of sequences into a matrix where each sequence is a row
  seqMatrix <- do.call(rbind, lapply(seqs, function(seq) strsplit(seq, "")[[1]]))
  
  # Find columns that do not contain gaps (assuming "-" as the gap symbol)
  validColumns <- !apply(seqMatrix, 2, function(column) any(column == "-"))
  
  # Extract these columns to create a cleaned matrix
  cleanedMatrix <- seqMatrix[, validColumns]
  
  # Convert cleaned matrix back to a list of sequences
  cleanedSeqs <- apply(cleanedMatrix, 1, paste, collapse = "")
  
  return(cleanedSeqs)
}


# -------------------------------------------------------------------------
# Count SNPs
countSNPsHelper <- function(fastaData, gapDeletion = TRUE) {
  refs <- fastaData$sequences
  refHeaders <- fastaData$headers
  
  # Optionally remove sites with missing data
  if (gapDeletion) {
    refs = deleteMissingDataSites(refs)
  }
  
  # Convert all cleaned sequences to a matrix of character vectors
  seqMatrix <- do.call(rbind, lapply(refs, function(seq) strsplit(seq, "")[[1]]))
  
  # Initialize the SNP matrix
  snpMatrix <- matrix(nrow = nrow(seqMatrix), ncol = nrow(seqMatrix), dimnames = list(refHeaders, refHeaders))
  
  # Loop over all pairs of sequences using matrix indices
  for (i in 1:nrow(seqMatrix)) {
    for (j in i:nrow(seqMatrix)) {
      seq1 <- seqMatrix[i, ]
      seq2 <- seqMatrix[j, ]
      
      # Calculate SNPs
      if (length(seq1) == length(seq2) && length(seq1) > 0) {
        snps <- sum(seq1 != seq2)
        snpMatrix[i, j] <- snpMatrix[j, i] <- snps
      } else {
        snpMatrix[i, j] <- snpMatrix[j, i] <- NA
      }
    }
  }
  
  return(snpMatrix)
}


# -------------------------------------------------------------------------
# p-distance
calcPDistance <- function(fastaData, gapDeletion = TRUE) {
  # Count the SNPs between each query and each reference sequence
  snpCounts <- countSNPsHelper(fastaData, gapDeletion = gapDeletion)
  
  # Extract query headers
  queryHeaders <- refHeaders <- fastaData$headers
  
  # Directly calculate reference sequence lengths
  refs <- fastaData$sequences

  # Optionally remove sites with missing data
  if (gapDeletion) {
    refs <- deleteMissingDataSites(refs)
  }
  
  refLengths <- sapply(refs, function(seq) length(unlist(strsplit(seq, ""))))
  
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

# -------------------------------------------------------------------------
# Jukes Cantor
# Function to calculate Jukes-Cantor genetic distance
calcJukesCantorDistance <- function(fastaData, gapDeletion = TRUE) {
  # Calculate p-distance for multiple queries
  p_dist <- calcPDistance(fastaData, gapDeletion = gapDeletion)  
  
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

# -------------------------------------------------------------------------
# helper function for Kimura2P and Tamura3P
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

# -------------------------------------------------------------------------

# Kimura 2 parameter
## Function to calculate Kimura 2-parameter genetic distance
calcKimura2pDistance <- function(fastaData, gapDeletion = TRUE) {
  
  # Load reference & query seq data and headers
  refs <- queries <- fastaData$sequences
  refHeaders <- queryHeaders <- fastaData$headers
  
  # Optionally remove sites with missing data
  if (gapDeletion) {
    refs <- queries <- deleteMissingDataSites(refs)
  }
  
  
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


# -------------------------------------------------------------------------

# calculate Tamura 3-parameter genetic distance
calcTamura3pDistance <- function(fastaData, gapDeletion = TRUE) {
  # Load reference & query seq data and headers
  refs <- queries <- fastaData$sequences
  refHeaders <- queryHeaders <- fastaData$headers
  
  # Optionally remove sites with missing data
  if (gapDeletion) {
    refs <- queries <- deleteMissingDataSites(refs)
  }
  
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
