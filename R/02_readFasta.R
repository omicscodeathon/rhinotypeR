# Function 1 
# Compare lengths among input sequences to pad short sequences with "-" until they are as long as the longest seq
compareLengths <- function(seqs) {
  # Calculate the maximum length from the list of sequences
  max_length <- max(sapply(seqs, nchar))
  
  # Adjust each sequence to have the same length as the longest sequence
  adjustedSeqs <- lapply(seqs, function(seq) {
    seq_length <- nchar(seq)
    if (seq_length < max_length) {
      # Pad shorter sequences with hyphens to reach the max length
      seq <- paste0(seq, paste0(rep("-", max_length - seq_length), collapse = ""))
    }
    seq
  })
  
  return(adjustedSeqs)
}



# Function 2
# Function to read sequences from a FASTA file and adjust their lengths
readFasta2 <- function(fastaFile) {
  # Read all lines from the FASTA file
  #lines <- readLines(fastaFile)
  
  alignment <- Biostrings::readDNAMultipleAlignment(fastaFile)
  
  # Initialize lists to store sequences and their headers
  # seqList <- list()
  # headerList <- c()
  
  seqList <- as.character(alignment) 
  headerList <- rownames(alignment)
  # 
  # # Temporary storage for the current sequence being read
  # currentSeq <- NULL
  # 
  # # Iterate through each line of the FASTA file
  # for (line in lines) {
  #   if (startsWith(line, ">")) {
  #     # If currentSeq is not NULL, it means we've finished reading a sequence
  #     if (!is.null(currentSeq)) {
  #       # Join all parts of the sequence into one
  #       fullSeq <- paste(currentSeq, collapse = "")
  #       seqList[[length(seqList) + 1]] <- fullSeq
  #     }
  #     # Reset currentSeq for the next sequence
  #     currentSeq <- c()
  #     # Add the header (without the ">" character) to headerList
  #     headerList <- c(headerList, substring(line, 2))
  #   } else {
  #     # If the line is not a header, it's part of the current sequence
  #     currentSeq <- c(currentSeq, toupper(line))
  #   }
  # }
  # 
  # # Add the last sequence to seqList if it exists
  # if (!is.null(currentSeq)) {
  #   fullSeq <- paste(currentSeq, collapse = "")
  #   seqList[[length(seqList) + 1]] <- fullSeq
  # }
  # 
  # Adjust all sequences to the length of the longest sequence
  seqList <- compareLengths(seqList)
  
  # Return a list containing the sequences and their corresponding headers
  return(list(sequences = seqList, headers = headerList))
}


# readFasta("inst/extdata/test.fasta")
# readFasta2("inst/extdata/input_aln.fasta")
# 
# 
# 
# # Read the DNA sequences from a FASTA file
# dna_sequences <- readDNAStringSet("./inst/extdata/test.fasta")
# # Translate the DNA sequences to amino acids
# amino_acid_sequences <- translate(dna_sequences)
# 
# seqList <- as.character(alignment) 
# headerList <- rownames(alignment)
