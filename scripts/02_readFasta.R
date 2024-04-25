
# Function 1 
# Compare lengths among input sequences to truncate very long sequences or pad short sequences with "-"
compareLengths <- function(seq, desired_length) {
  seq_length <- nchar(seq)
  # If the sequence is longer than 430 characters, truncate it to 430 characters
  if (seq_length > desired_length){
    seq <- substr(seq, 1,desired_length)
    
    # If the sequence is shorter than 430 characters, pad it with hyphens
  }else if (seq_length < desired_length) {
    seq <- paste0(seq, paste0(rep("-", desired_length - seq_length), collapse = ""))
  }
  return(seq)
}




# Function 2
# Function to read sequences from a FASTA file and adjust their lengths
readFasta <- function(fastaFile, desiredLength = 430) {
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
      if (!is.null(currentSeq)) {
        # Join all parts of the sequence into one and adjust its length
        fullSeq <- paste(currentSeq, collapse = "")
        adjustedSeq <- compareLengths(fullSeq, desiredLength)
        seqList[[length(seqList) + 1]] <- adjustedSeq
      }
      # Reset currentSeq for the next sequence
      currentSeq <- c()
      # Add the header (without the ">" character) to headerList
      headerList <- c(headerList, substring(line, 2))
    } else {
      # If the line is not a header, it's part of the current sequence
      currentSeq <- c(currentSeq, toupper(line))
    }
  }
  
  # After the loop, add the last sequence to seqList if it exists
  if (!is.null(currentSeq)) {
    fullSeq <- paste(currentSeq, collapse = "")
    adjustedSeq <- compareLengths(fullSeq, desiredLength)
    seqList[[length(seqList) + 1]] <- adjustedSeq
  }
  
  # Return a list containing the sequences and their corresponding headers
  return(list(sequences = seqList, headers = headerList))
}


# Example usage
# readFasta("./data/RVAPrototypeAligned.fasta")
# readFasta("./data/tmp.fasta")
