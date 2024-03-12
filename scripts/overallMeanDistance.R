


# Function to calculate overall mean distance of a multiple sequence alignment
overallMeanDistance <- function(fastaFile, model ='p-distance') {
  # First, read the sequences from the FASTA file
  fastaData <- readFasta(fastaFile)
  sequences <- fastaData$sequences
  
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
  overall_mean_distance <- total_distance / num_comparisons
  
  return(overall_mean_distance)
}


# Example usage
overallMeanDistance("./data/RVAPrototypeAligned.fasta")  




