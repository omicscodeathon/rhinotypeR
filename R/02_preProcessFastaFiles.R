# Function 1 
# Compare lengths among input sequences to pad short sequences with "-" until they are as long as the longest seq
compareLengths <- function(seqs) {
  # Calculate the maximum length from the list of sequences
  max_length <- max(nchar(seqs))
  
  # Adjust each sequence to have the same length as the longest sequence
  adjustedSeqs <- vapply(seqs, function(seq) {
    seq_length <- nchar(seq)
    if (seq_length < max_length) {
      # Pad shorter sequences with hyphens to reach the max length
      seq <- paste0(seq, paste0(rep("-", max_length - seq_length), collapse = ""))
    }
    seq
  }, character(1), USE.NAMES = TRUE)
  
  return(adjustedSeqs)
}


# Read the DNA sequences from a FASTA file (Biostrings package)
# readDNAStringSetObject <- Biostrings::readDNAStringSet(fastaFile)


# Function 2
# Function to read sequences from a FASTA file and adjust their lengths
preProcessFastaStringSet <- function(readDNAStringSetObject) {
  
  # Convert the alignment to character sequences and extract headers
  seqList <- as.character(readDNAStringSetObject) 
  headerList <- names(seqList)

  # Adjust all sequences to the length of the longest sequence
  seqList <- compareLengths(seqList)
  
  # Return a list containing the sequences and their corresponding headers
  return(list(sequences = seqList, headers = headerList))
}
