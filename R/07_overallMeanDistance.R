
# Function to calculate overall mean distance of a multiple sequence alignment

source("./scripts/04_genetic_distances.R")

 # 1. p-distance

overallPDistance <- function(fastaData, gapDeletion=TRUE) {
   # extract seq data from the fastaData (a product of readFasta)
  sequences <- fastaData$sequences
  
  # Optionally remove sites with missing data
  if (gapDeletion) {
    sequences <- deleteMissingDataSites(sequences)
  }
  
  num_sequences <- length(sequences)
  total_distance <- 0
  num_comparisons <- 0
  
  for (i in 1:(num_sequences - 1)) {
    for (j in (i + 1):num_sequences) {
      # Compute distance between sequence i and j
      seq_i_chars <- strsplit(sequences[[i]], "")[[1]]
      seq_j_chars <- strsplit(sequences[[j]], "")[[1]]
      distance <- sum(seq_i_chars != seq_j_chars) / length(seq_i_chars)
      total_distance <- total_distance + distance
      num_comparisons <- num_comparisons + 1
    }
  }
  
  # Calculate overall mean distance
  overallPDistance <- total_distance / num_comparisons
  
  return(overallPDistance)
}

# -------------------------------------------------------------------------

# 2. Jukes-Cantor model
overallJCDistance <- function(fastaData, gapDeletion=TRUE) {
  sequences <- fastaData$sequences
  num_sequences <- length(sequences)
  total_jc_distance <- 0
  num_comparisons <- 0
  
  # Optionally remove sites with missing data
  if (gapDeletion) {
    sequences <- deleteMissingDataSites(sequences)
  }
  
  for (i in 1:(num_sequences - 1)) {
    for (j in (i + 1):num_sequences) {
      # Compute distance between sequence i and j using p-distance
      seq_i_chars <- strsplit(sequences[[i]], "")[[1]]
      seq_j_chars <- strsplit(sequences[[j]], "")[[1]]
      p_distance <- sum(seq_i_chars != seq_j_chars) / length(seq_i_chars)
      
      # Apply the Jukes-Cantor correction if p_distance < 3/4, else consider the distance as infinite
      if (p_distance < 0.75) {
        jc_distance <- -3/4 * log(1 - 4/3 * p_distance)
        total_jc_distance <- total_jc_distance + jc_distance
        num_comparisons <- num_comparisons + 1
      }
    }
  }
  
  # Calculate overall mean Jukes-Cantor distance
  if (num_comparisons > 0) {
    overall_mean_jc_distance <- total_jc_distance / num_comparisons
  } else {
    overall_mean_jc_distance <- NA  # Return NA if no valid comparisons were made
  }
  
  return(overall_mean_jc_distance)
}

# -------------------------------------------------------------------------


# Kimura 2 parameter
overallK2PDistance <- function(fastaData, gapDeletion=TRUE) {
  sequences <- fastaData$sequences
  
  num_sequences <- length(sequences)
  total_distance <- 0
  num_comparisons <- 0

  # Optionally remove sites with missing data
  if (gapDeletion) {
    sequences <- deleteMissingDataSites(sequences)
  }
  
  
  for (i in 1:(num_sequences - 1)) {
    for (j in (i + 1):num_sequences) {
      seq_i_chars <- strsplit(sequences[[i]], "")[[1]]
      seq_j_chars <- strsplit(sequences[[j]], "")[[1]]
      
      # Determine transitions and transversions
      transitions <- 0
      transversions <- 0
      for (k in 1:length(seq_i_chars)) {
        if (seq_i_chars[k] != seq_j_chars[k]) {
          # Purines: A, G; Pyrimidines: C, T
          if (seq_i_chars[k] %in% c("A", "G") && seq_j_chars[k] %in% c("A", "G") ||
              seq_i_chars[k] %in% c("C", "T") && seq_j_chars[k] %in% c("C", "T")) {
            transitions <- transitions + 1
          } else {
            transversions <- transversions + 1
          }
        }
      }
      
      P <- transitions / length(seq_i_chars)
      Q <- transversions / length(seq_i_chars)
      
      # Apply Kimura 2-parameter model if valid
      if ((1 - 2*P - Q) > 0 && (1 - 2*Q) > 0) {
        k2p_distance <- 0.5 * log(1 / (1 - 2*P - Q)) + 0.25 * log(1 / (1 - 2*Q))
        total_distance <- total_distance + k2p_distance
      } else {
        # cases where K2P model is not valid
      }
      
      num_comparisons <- num_comparisons + 1
    }
  }
  
  if (num_comparisons > 0) {
    overall_mean_distance <- total_distance / num_comparisons
  } else {
    overall_mean_distance <- NA # No valid comparisons or all were skipped due to invalid K2P model conditions
  }
  
  return(overall_mean_distance)
}

# -------------------------------------------------------------------------

# Tamura 3 parameter

overallT3PDistance <- function(fastaData, gapDeletion=TRUE) {
  sequences <- fastaData$sequences
  
  num_sequences <- length(sequences)
  total_distance <- 0
  num_comparisons <- 0
  
  # Optionally remove sites with missing data
  if (gapDeletion) {
    sequences <- deleteMissingDataSites(sequences)
  }
  
  for (i in 1:(num_sequences - 1)) {
    for (j in (i + 1):num_sequences) {
      seq_i_chars <- strsplit(sequences[[i]], "")[[1]]
      seq_j_chars <- strsplit(sequences[[j]], "")[[1]]
      
      # Use existing function to calculate P and Q
      transitionTransversionProportions <- calcTransitionTransversions(seq_i_chars, seq_j_chars)
      P <- transitionTransversionProportions$P
      Q <- transitionTransversionProportions$Q
      
      gc_content <- sum(seq_i_chars %in% c("G", "C")) / length(seq_i_chars)
      
      # Apply Tamura 3-parameter model if valid
      if ((1 - 2*P - Q) > 0 && (1 - 2*Q) > 0) {
        G <- 1 / (1 + (1 - 2*Q) * (gc_content / (1 - gc_content)))
        t3p_distance <- (0.5 * log(1 / (1 - 2*P - Q))) + (G * 0.5 * log(1 / (1 - 2*Q)))
        total_distance <- total_distance + t3p_distance
      } else {
        # Optionally, handle cases where T3P model is not valid
        # For this example, we'll skip these cases in the average calculation
      }
      
      num_comparisons <- num_comparisons + 1
    }
  }
  
  if (num_comparisons > 0) {
    overall_mean_distance <- total_distance / num_comparisons
  } else {
    overall_mean_distance <- NA # skipped due to invalid T3P model conditions
  }
  
  return(overall_mean_distance)
}

# -------------------------------------------------------------------------
# Function to bring all evo models together
overallMeanDistance<- function(fastaData, model, gapDeletion=TRUE) {
  
  # Determine which model to use based on user input
  if (model == "p-distance") {
    result <- overallPDistance(fastaData, gapDeletion = gapDeletion)
  } else if (model == "JC") {
    result <- overallJCDistance(fastaData, gapDeletion = gapDeletion)
  } else if (model == "Kimura2p") {
    result <- overallK2PDistance(fastaData, gapDeletion = gapDeletion)
  } else if (model == "Tamura3p") {
    result <- overallT3PDistance(fastaData, gapDeletion = gapDeletion)
  } else {
    stop("Unknown model specified. Choose from 'p-distance', 'JC', 'Kimura2p', or 'Tamura3p' ")
  }
  
  # Return the result of the chosen model
  return(result)
}
