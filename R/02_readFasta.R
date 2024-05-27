# Function 1 
# Compare lengths among input sequences to pad short sequences with "-" until they are as long as the longest seq
compareLengths <- function(seqs) {
  # Calculate the maximum length from the list of sequences
  max_length <- max(nchar(seqs))
  
  # Adjust each sequence to have the same length as the longest sequence
  adjustedSeqs <- sapply(seqs, function(seq) {
    seq_length <- nchar(seq)
    if (seq_length < max_length) {
      # Pad shorter sequences with hyphens to reach the max length
      seq <- paste0(seq, paste0(rep("-", max_length - seq_length), collapse = ""))
    }
    seq
  }, USE.NAMES = TRUE)
  
  return(adjustedSeqs)
}



# Function 2
# Function to read sequences from a FASTA file and adjust their lengths
readFasta <- function(fastaFile) {
  # Read the DNA sequences from a FASTA file
  alignment <- Biostrings::readDNAMultipleAlignment(fastaFile)
  
  # Convert the alignment to character sequences and extract headers
  seqList <- as.character(alignment) 
  headerList <- rownames(alignment)
  
  # Combine headers and sequences into a named character vector
  names(seqList) <- headerList
  
  # Adjust all sequences to the length of the longest sequence
  seqList <- compareLengths(seqList)
  
  # Return a list containing the sequences and their corresponding headers
  return(list(sequences = seqList, headers = headerList))
}
