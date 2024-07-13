source("R/02_readFasta.R")


readAA <- function(fastaFile) {
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
        # Join all parts of the sequence into one
        fullSeq <- paste(currentSeq, collapse = "")
        seqList[[length(seqList) + 1]] <- fullSeq
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

  # Add the last sequence to seqList if it exists
  if (!is.null(currentSeq)) {
    fullSeq <- paste(currentSeq, collapse = "")
    seqList[[length(seqList) + 1]] <- fullSeq
  }

  # Adjust all sequences to the length of the longest sequence
  seqList <- compareLengths(seqList)
  
  # Return a list containing the sequences and their corresponding headers
  return(list(sequences = seqList, headers = headerList))
}



plotAA <- function(fastaFile, showLegend = FALSE) {
  fastaData <- readAA(fastaFile = fastaFile)  # Read file
  
  sequences <- fastaData$sequences
  seqNames <- fastaData$headers
  proteinLength <- max(vapply(sequences, nchar, integer(1)))
  
  compareSequences <- function(seqA, seqB) {
    seqAChars <- strsplit(seqA, "")[[1]]
    seqBChars <- strsplit(seqB, "")[[1]]
    differences <- which(seqAChars != seqBChars)
    subsType <- seqBChars[differences]
    data.frame(position = differences, subsType = subsType)
  }
  
  colorMap <- c(
    R = "red", H = "red", K = "red",      # Positively charged amino acid 
    D = "blue", E = "blue",               # Negatively charged amino acid 
    S = "green", T = "green", N = "green", Q = "green", # Polar amino acid 
    A = "yellow", V = "yellow", I = "yellow", L = "yellow", M = "yellow", F = "yellow", 
    W = "yellow", P = "yellow", G = "yellow", Y = "yellow", C = "yellow"  # Nonpolar amino acid 
  )
  
  diffList <- list()
  for (i in 2:length(sequences)) {
    diffList[[i - 1]] <- compareSequences(sequences[[1]], sequences[[i]])
    diffList[[i - 1]]$color <- colorMap[diffList[[i - 1]]$subsType]
    isOther <- is.na(diffList[[i - 1]]$color)
    diffList[[i - 1]]$color[isOther] <- "gray"
  }
  
  # Adjust left margin to ensure y-axis labels are not truncated
  oldPar <- par(mar = c(5, 8, 4, 2) + 0.1, xpd = TRUE) # Increase the left margin and allow plotting outside plot
  plot(NULL, xlim = c(1, proteinLength), ylim = c(0.5, length(sequences)), type = 'n',
       xlab = paste0("Protein Position of ", seqNames[length(seqNames)], ", acting as reference"),
       ylab = "", 
       yaxt = 'n')
  axis(2, at = seq_along(sequences), labels = seqNames, las = 2, cex.axis = 0.8)
  
  # Plot small vertical bars for each difference using mapped colors
  for (i in seq_along(diffList)) {
    yPosStart <- rep(i, nrow(diffList[[i]]))
    yPosEnd <- yPosStart
    for (j in seq_len(nrow(diffList[[i]]))) {
      segments(x0 = diffList[[i]]$position[j], y0 = yPosStart[j] - 0.25, 
               x1 = diffList[[i]]$position[j], y1 = yPosEnd[j] + 0.25, col = diffList[[i]]$color[j])
    }
  }
  if (showLegend){ # Add a semi-transparent legend in the top-left corner
  legend("topleft", inset=c(0.78, 0),
         legend=c("+ve charged", "-ve charged", "Polar", "Non-polar", "Other"),
         fill=c("red", "blue", "green", "yellow", "gray"), cex=0.45, bty="n", 
         box.col="gray", bg=adjustcolor("white", alpha.f=0.7))
  
  par(oldPar) # Reset to old graphical parameters after plotting
  }
}
