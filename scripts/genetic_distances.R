

# Base R functions


# Function 1 
# Function to read sequences from a FASTA file
readFasta <- function(fastaFile) {
  # Read all lines from the FASTA file
  lines <- readLines(fastaFile)
  
  # Initialize lists to store sequences and their headers
  seqList <- list()
  headerList <- c()
  
  # Temporary storage for the current sequence being read
  currentSeq <- NULL
  
  # Iterate through each line of the FASTA file
  for (line in lines) {
    if (startsWith(line, ">")) {
      # If currentSeq is not NULL, it means we've finished reading a sequence
      # Add it to seqList
      if (!is.null(currentSeq)) {
        seqList[[length(seqList) + 1]] <- paste(currentSeq, collapse = "")
      }
      # Reset currentSeq for the next sequence
      currentSeq <- c()
      # Add the header (without the ">" character) to headerList
      headerList <- c(headerList, substring(line, 2))
    } else {
      # If the line is not a header, it's part of the current sequence
      # Convert it to uppercase and add it to currentSeq
      currentSeq <- c(currentSeq, toupper(line))
    }
  }
  
  # After the loop, add the last sequence to seqList if it exists
  if (!is.null(currentSeq)) {
    seqList[[length(seqList) + 1]] <- paste(currentSeq, collapse = "")
  }
  
  # Return a list containing the sequences and their corresponding headers
  return(list(sequences = seqList, headers = headerList))
}




# Example usage
readFasta("./data/RVAPrototypeAligned.fasta")
readFasta("./data/tmp_ref_b99.fasta")

# Function 2
# Compare sequence lengths
# Modified function to compare lengths and calculate the difference between multiple references and one query sequence
# Adjusted function to compare lengths and calculate the difference between references and queries using sequence headers for naming
compareLengths <- function(pathToRef, pathToQuery) {
  # Read multiple reference sequences and their headers
  fastaRef <- readFasta(pathToRef)
  refs <- fastaRef$sequences
  refHeaders <- fastaRef$headers
  
  # Read query sequences and their headers
  fastaQuery <- readFasta(pathToQuery)
  queries <- fastaQuery$sequences
  queryHeaders <- fastaQuery$headers
  
  # Assuming the first query for simplicity in this example; adjust as needed for multiple queries
  query <- queries[[1]]
  queryHeader <- queryHeaders[[1]]
  
  # Initialize a matrix to store the differences in lengths
  # Single row for the query and columns for each reference, using headers for naming
  lengthMatrix <- matrix(nrow = 1, ncol = length(refs),
                         dimnames = list(c(queryHeader), refHeaders))
  
  # Calculate the query length
  queryLength <- length(unlist(strsplit(query, "")))
  
  # Populate the matrix with the difference in length between the query and each reference sequence
  for (i in seq_along(refs)) {
    refLength <- length(unlist(strsplit(refs[[i]], "")))
    
    # Store the difference in length
    lengthMatrix[1, i] <- refLength - queryLength
  }
  
  # Return the matrix
  return(lengthMatrix)
}


# Example usage
compareLengths(pathToRef = "./data/RVAPrototypeAligned.fasta", pathToQuery="./data/tmp_query.fasta")

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
countSNPs <- function(pathToRef, pathToQuery) {
  # Read refs and query sequences along with their headers
  fastaRef <- readFasta(pathToRef)
  refs <- fastaRef$sequences
  refHeaders <- fastaRef$headers
  
  fastaQuery <- readFasta(pathToQuery)
  query <- fastaQuery$sequences[[1]] # Assuming only one query sequence
  queryHeader <- fastaQuery$headers[[1]] # Assuming only one query header
  
  # Prepare the matrix to store SNP counts, using the query header for row name
  snpMatrix <- matrix(nrow = 1, ncol = length(refs), dimnames = list(c(queryHeader), refHeaders))
  
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
    snpMatrix[1, i] <- snps
  }
  
  # Output
  return(snpMatrix)
}



# Example usage
countSNPs(pathToRef = "./data/RVBPrototypeAligned.fasta", pathToQuery="./data/tmp_query.fasta")



# Function 5
# p-distance
calcPDistance <- function(pathToRef, pathToQuery) {
  # Count the SNPs between the query and each reference sequence
  snpCounts <- countSNPs(pathToRef, pathToQuery)
  
  # Read query sequence to obtain the query header
  fastaQuery <- readFasta(pathToQuery)
  queryHeader <- fastaQuery$headers[1]  # Assuming only one query sequence
  
  # Directly calculate reference sequence lengths
  fastaRef <- readFasta(pathToRef)
  refHeaders <- fastaRef$headers
  refLengths <- sapply(fastaRef$sequences, function(seq) length(unlist(strsplit(seq, ""))))
  
  # Prepare a matrix for p-distances with appropriate dimensions and names
  pDistancesMatrix <- matrix(nrow = 1, ncol = length(refLengths), dimnames = list(queryHeader, refHeaders))
  
  # Calculate p-distance for each reference sequence
  for (i in seq_along(refLengths)) {
    if (!is.na(snpCounts[1, i])) {
      pDistancesMatrix[1, i] <- snpCounts[1, i] / refLengths[i]
    } else {
      pDistancesMatrix[1, i] <- NA  # Assign NA if SNP count was NA
    }
  }
  
  return(pDistancesMatrix)
}


# Example usage
calcPDistance(pathToRef = "./data/RVBPrototypeAligned.fasta", pathToQuery = "./data/tmp_query.fasta")


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
calcJukesCantorDistance(pathToRef = "./data/RVBPrototypeAligned.fasta", pathToQuery = "./data/tmp_query.fasta")


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
  # Read multiple reference sequences and one query sequence
  refsData <- readFasta(pathToRef)
  refs <- refsData$sequences
  refHeaders <- refsData$headers
  
  queryData <- readFasta(pathToQuery)
  query <- queryData$sequences[[1]]  # Assuming only one query sequence
  queryHeader <- queryData$headers[[1]]  # Assuming only one query header
  
  # Initialize a matrix to store the distances
  k2pMatrix <- matrix(nrow = 1, ncol = length(refs), 
                      dimnames = list(c(queryHeader), refHeaders))
  
  # Convert the query sequence to character vector
  queryChars <- strsplit(query, split = "")[[1]]
  
  # Iterate over each reference sequence
  for (i in seq_along(refs)) {
    refChars <- strsplit(refs[[i]], split = "")[[1]]
    
    # Calculate proportions of transitions and transversions
    tt <- calcTransitionTransversions(refChars, queryChars)
    P <- tt$P
    Q <- tt$Q
    
    # Calculate Kimura 2-parameter genetic distance
    K2P_distance <- -0.5 * log((1 - 2*P - Q) * sqrt(1 - 2*Q))
    
    # Store the distance
    k2pMatrix[1, i] <- K2P_distance
  }
  
  # Return the matrix of Kimura 2-parameter genetic distances
  return(k2pMatrix)
}


# Example usage
calcKimura2pDistance(pathToRef = "./data/RVBPrototypeAligned.fasta", pathToQuery = "./data/tmp_query.fasta")


# calculate Tamura 3-parameter genetic distance
calcTamura3pDistance <- function(pathToRef, pathToQuery) {
  # Read multiple reference sequences and one query sequence
  refsData <- readFasta(pathToRef)
  refs <- refsData$sequences
  refHeaders <- refsData$headers
  
  queryData <- readFasta(pathToQuery)
  query <- queryData$sequences[[1]]  # Assuming only one query sequence
  queryHeader <- queryData$headers[[1]]  # Assuming only one query header
  
  # Initialize a matrix to store the distances
  tamuraMatrix <- matrix(nrow = 1, ncol = length(refs), 
                         dimnames = list(c(queryHeader), refHeaders))
  
  # Convert the query sequence to character vector
  queryChars <- strsplit(query, split = "")[[1]]
  
  # Iterate over each reference sequence
  for (i in seq_along(refs)) {
    refChars <- strsplit(refs[[i]], split = "")[[1]]
    
    # Calculate base frequencies for both sequences
    gcContentRef <- sum(refChars %in% c("G", "C")) / length(refChars)
    gcContentQuery <- sum(queryChars %in% c("G", "C")) / length(queryChars)
    
    # Calculate transitions and transversions
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
    
    # Store the distance
    tamuraMatrix[1, i] <- distance
  }
  
  # Return the matrix of Tamura 3-parameter genetic distances
  return(tamuraMatrix)
}


# Example usage
calcTamura3pDistance(pathToRef = "./data/RVBPrototypeAligned.fasta", pathToQuery = "./data/tmp_query.fasta")


# Tamura-Nei model  (complex!!!)

# General Time Reversible (GTR) Model (complex!!!)

# Maximum Composite Likelihood method? (possibly too complex for regular scripting!!!)





allPrototypeDistances <- function(pathToRef, pathToQuery, model = "p-distance") {
  # Determine which model to use based on user input
  if (model == "p-distance") {
    result <- calcPDistance(pathToRef, pathToQuery)
  } else if (model == "JC") {
    result <- calcJukesCantorDistance(pathToRef, pathToQuery)
  } else if (model == "Kimura2p") {
    result <- calcKimura2pDistance(pathToRef, pathToQuery)
  } else if (model == "Tamura3p") {
    result <- calcTamura3pDistance(pathToRef, pathToQuery)
  } else {
    stop("Unknown model specified. Choose from 'p-distance', 'JC', 'Kimura2p', or 'Tamura3p'.")
  }
  
  # Return the result of the chosen model
  return(result)
}



# Example usage
allPrototypeDistances("./data/RVBPrototypeAligned.fasta", "./data/tmp_query.fasta", "p-distance")
allPrototypeDistances("./data/RVBPrototypeAligned.fasta", "./data/tmp_query.fasta", "JC")
allPrototypeDistances("./data/RVBPrototypeAligned.fasta", "./data/tmp_query.fasta", "Kimura2p")
allPrototypeDistances("./data/RVBPrototypeAligned.fasta", "./data/tmp_query.fasta", "Tamura3p")




