

# Optionally delete sites with at least one missing data for all sequences 
deleteMissingDataSites <- function(seqs) {
  # Convert list of sequences into a matrix where each sequence is a row
  seqMatrix <- do.call(rbind, lapply(seqs, function(seq) strsplit(seq, "")[[1]]))
  
  # Find columns that do not contain gaps (assuming "-" as the gap symbol)
  validColumns <- !apply(seqMatrix, 2, function(column) any(column == "-"))
  
  # Extract these columns to create a cleaned matrix, drop = FALSE to avoid reduction to a vector
  cleanedMatrix <- seqMatrix[, validColumns, drop = FALSE]
  
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
    refs <- deleteMissingDataSites(refs)
  }
  
  # Convert all cleaned sequences to a matrix of character vectors
  seqMatrix <- do.call(rbind, lapply(refs, function(seq) strsplit(seq, "")[[1]]))
  
  # Function to calculate SNPs between two sequences
  calcSNPs <- function(seq1, seq2) {
    if (length(seq1) == length(seq2) && length(seq1) > 0) {
      return(sum(seq1 != seq2))
    } else {
      return(NA)
    }
  }
  
  # Apply the SNP calculation function over all pairs of sequences
  snpMatrix <- outer(
    seq_len(nrow(seqMatrix)), 
    seq_len(nrow(seqMatrix)), 
    Vectorize(function(i, j) calcSNPs(seqMatrix[i, ], seqMatrix[j, ]))
  )
  
  # Set the dimnames of the matrix
  dimnames(snpMatrix) <- list(refHeaders, refHeaders)
  
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
  
  refLengths <- vapply(refs, function(seq) length(unlist(strsplit(seq, ""))), integer(1))
  
  # Vectorized calculation of p-distances for each query-reference pair
  pDistancesMatrix <- outer(
    seq_len(nrow(snpCounts)),
    seq_along(refLengths),
    Vectorize(function(q, i) {
      if (!is.na(snpCounts[q, i])) {
        return(snpCounts[q, i] / refLengths[i])
      } else {
        return(NA)
      }
    })
  )
  
  # Set the dimension names for the resulting matrix
  dimnames(pDistancesMatrix) <- list(queryHeaders, refHeaders)
  
  return(pDistancesMatrix)
}



# -------------------------------------------------------------------------
# Jukes Cantor
# Function to calculate Jukes-Cantor genetic distance
calcJukesCantorDistance <- function(fastaData, gapDeletion = TRUE) {
  # Calculate p-distance for multiple queries
  p_dist <- calcPDistance(fastaData, gapDeletion = gapDeletion)  
  
  # Initialize a matrix to store Jukes-Cantor distances
  jc_dist <- matrix(nrow = nrow(p_dist), 
                    ncol = ncol(p_dist), 
                    dimnames = dimnames(p_dist))
  
  # Apply the Jukes-Cantor formula to each element in the p-distance matrix
  jc_dist <- -3/4 * log(1 - 4/3 * p_dist)
  
  # Handling cases where p_dist >= 0.75, setting JC distance to Inf 
  # JC assumes that all nt substitutions are equally probable and independent, which might not always hold true in real data. 
      #The handling of cases where p_dist >= 0.75 with Inf is to indicate that the JC model might not be valid for these 
        #high levels of divergence due to the assumption of the model being violated. 
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
  
  # Load reference & query sequence data and headers
  refs <- queries <- fastaData$sequences
  refHeaders <- queryHeaders <- fastaData$headers
  
  # Optionally remove sites with missing data
  if (gapDeletion) {
    refs <- queries <- deleteMissingDataSites(refs)
  }
  
  # Define a function to calculate the K2P distance for a pair of sequences
  k2p_distance_function <- function(query, ref) {
    queryChars <- strsplit(query, split = "")[[1]]
    refChars <- strsplit(ref, split = "")[[1]]
    
    # Ensure lengths match to avoid errors in distance calculation
    if (length(queryChars) == length(refChars)) {
      # Calculate proportions of transitions and transversions
      tt <- calcTransitionTransversions(refChars, queryChars)
      P <- tt$P
      Q <- tt$Q
      
      # Calculate Kimura 2-parameter genetic distance
      if ((1 - 2 * P - Q) > 0 && (1 - 2 * Q) > 0) {
        return(-0.5 * log((1 - 2 * P - Q) * sqrt(1 - 2 * Q)))
      } else {
        return(NA)  # Return NA if the calculation is not possible
      }
    } else {
      return(NA)  # Return NA if sequence lengths do not match
    }
  }
  
  # Vectorized calculation of K2P distances
  k2pMatrix <- outer(queries, refs, Vectorize(k2p_distance_function))
  
  # Set the dimension names for the resulting matrix
  dimnames(k2pMatrix) <- list(queryHeaders, refHeaders)
  
  # Return the matrix of Kimura 2-parameter genetic distances
  return(k2pMatrix)
}



# -------------------------------------------------------------------------

# calculate Tamura 3-parameter genetic distance
calcTamura3pDistance <- function(fastaData, gapDeletion = TRUE) {
  # Load reference & query sequence data and headers
  refs <- queries <- fastaData$sequences
  refHeaders <- queryHeaders <- fastaData$headers
  
  # Optionally remove sites with missing data
  if (gapDeletion) {
    refs <- queries <- deleteMissingDataSites(refs)
  }
  
  # Define a function to calculate the Tamura 3-parameter distance for a pair of sequences
  tamura3p_distance_function <- function(query, ref) {
    queryChars <- strsplit(query, split = "")[[1]]
    refChars <- strsplit(ref, split = "")[[1]]
    
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
      if ((1 - P / C - Q) > 0 && (1 - 2 * Q) > 0) {
        return(-C * log(1 - P / C - Q) - 0.5 * (1 - C) * log(1 - 2 * Q))
      } else {
        return(NA)  # Assign NA if distance cannot be calculated
      }
    } else {
      return(NA)  # Assign NA if sequence lengths do not match
    }
  }
  
  # Vectorized calculation of Tamura 3-parameter distances
  tamuraMatrix <- outer(queries, refs, Vectorize(tamura3p_distance_function))
  
  # Set the dimension names for the resulting matrix
  dimnames(tamuraMatrix) <- list(queryHeaders, refHeaders)
  
  # Return the matrix of Tamura 3-parameter genetic distances
  return(tamuraMatrix)
}

