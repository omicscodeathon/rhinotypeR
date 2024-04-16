
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
#readFasta("./data/RVAPrototypeAligned.fasta")
#readFasta("./data/tmp_ref_b99.fasta")
