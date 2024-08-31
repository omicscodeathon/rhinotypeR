
# Function to calculate overall mean distance of a multiple sequence alignment

# Consolidated Helper function to reduce code repetition
calculateOverallDistancesHelper <- function(sequences, gapDeletion, distanceFunc) {
  # Optionally remove sites with missing data
  if (gapDeletion) {
    sequences <- deleteMissingDataSites(sequences)
  }
  
  # Calculate pairwise distances using the provided distance function
  pairwise_distances <- combn(sequences, 2, function(x) {
    seq_i_chars <- strsplit(x[1], "")[[1]]
    seq_j_chars <- strsplit(x[2], "")[[1]]
    
    distanceFunc(seq_i_chars, seq_j_chars)
  }, simplify = TRUE)
  
  # Calculate the overall mean distance
  overall_mean_distance <- mean(pairwise_distances, na.rm = TRUE)
  
  return(overall_mean_distance)
}



 # 1. p-distance

overallPDistance <- function(fastaData, gapDeletion = TRUE) {
  sequences <- fastaData$sequences
  
  distanceFunc <- function(seq_i_chars, seq_j_chars) {
    sum(seq_i_chars != seq_j_chars) / length(seq_i_chars)
  }
  
  calculateOverallDistancesHelper(sequences, gapDeletion, distanceFunc)
}



# -------------------------------------------------------------------------

# 2. Jukes-Cantor model
overallJCDistance <- function(fastaData, gapDeletion = TRUE) {
  sequences <- fastaData$sequences
  
  distanceFunc <- function(seq_i_chars, seq_j_chars) {
    p_distance <- sum(seq_i_chars != seq_j_chars) / length(seq_i_chars)
    if (p_distance < 0.75) {
      -3/4 * log(1 - 4/3 * p_distance)
    } else {
      NA
    }
  }
  
  calculateOverallDistancesHelper(sequences, gapDeletion, distanceFunc)
}



# -------------------------------------------------------------------------

# Kimura 2 parameter

overallK2PDistance <- function(fastaData, gapDeletion = TRUE) {
  sequences <- fastaData$sequences
  
  # Define the distance function for the Kimura 2-parameter model
  distanceFunc <- function(seq_i_chars, seq_j_chars) {
    transitions <- sum((seq_i_chars == "A" & seq_j_chars == "G") |
                         (seq_i_chars == "G" & seq_j_chars == "A") |
                         (seq_i_chars == "C" & seq_j_chars == "T") |
                         (seq_i_chars == "T" & seq_j_chars == "C"))
    
    transversions <- sum((seq_i_chars == "A" & seq_j_chars %in% c("C", "T")) |
                           (seq_i_chars == "G" & seq_j_chars %in% c("C", "T")) |
                           (seq_i_chars == "C" & seq_j_chars %in% c("A", "G")) |
                           (seq_i_chars == "T" & seq_j_chars %in% c("A", "G")))
    
    P <- transitions / length(seq_i_chars)
    Q <- transversions / length(seq_i_chars)
    
    # Apply the Kimura 2-parameter formula, ensuring the logarithm arguments are valid
    if ((1 - 2 * P - Q) > 0 && (1 - 2 * Q) > 0) {
      0.5 * log(1 / (1 - 2 * P - Q)) + 0.25 * log(1 / (1 - 2 * Q))
    } else {
      NA
    }
  }
  
  # Use the helper function to calculate the overall mean distance
  overall_mean_distance <- calculateOverallDistancesHelper(sequences, gapDeletion, distanceFunc)
  
  return(overall_mean_distance)
}


# -------------------------------------------------------------------------

# Tamura 3 parameter

overallT3PDistance <- function(fastaData, gapDeletion = TRUE) {
  sequences <- fastaData$sequences
  
  distanceFunc <- function(seq_i_chars, seq_j_chars) {
    tt_results <- calcTransitionTransversions(seq_i_chars, seq_j_chars)
    P <- tt_results$P
    Q <- tt_results$Q
    
    gc_content <- sum(seq_i_chars %in% c("G", "C")) / length(seq_i_chars)
    C <- gc_content + gc_content - 2 * gc_content * gc_content
    
    if ((1 - P / C - Q) > 0 && (1 - 2 * Q) > 0) {
      G <- 1 / (1 + (1 - 2 * Q) * (gc_content / (1 - gc_content)))
      -C * log(1 - P / C - Q) - 0.5 * (1 - C) * log(1 - 2 * Q)
    } else {
      NA
    }
  }
  
  calculateOverallDistancesHelper(sequences, gapDeletion, distanceFunc)
}


# -------------------------------------------------------------------------
# Function to bring all evo models together
overallMeanDistance <- function(fastaData, model = 'p-distance', gapDeletion = TRUE) {
  
  # Preprocess fasta data
  fastaData <- preProcessFastaStringSet(fastaData)
  
  # Map model names to their corresponding functions
  functionMap <- list(
    "p-distance" = overallPDistance,
    "JC" = overallJCDistance,
    "Kimura2p" = overallK2PDistance,
    "Tamura3p" = overallT3PDistance
  )
  
  # Apply the appropriate function based on the model
  result <- applyModelFunction(fastaData, model, gapDeletion, functionMap)
  
  return(result)
}

