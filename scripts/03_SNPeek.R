

source("./scripts/02_readFasta.R")

SNPeek <- function(fastaData, showLegend = FALSE) {
  sequences <- fastaData$sequences
  seqNames <- fastaData$headers
  genomeLength <- max(sapply(sequences, nchar))
  
  # Enhanced function to also return the substitution type
  compareSequences <- function(seqA, seqB) {
    seqAChars <- strsplit(seqA, "")[[1]]
    seqBChars <- strsplit(seqB, "")[[1]]
    differences <- which(seqAChars != seqBChars)
    subsType <- seqBChars[differences]
    data.frame(position = differences, subsType = subsType)
  }
  
  # Define a color map for substitutions
  colorMap <- c(A = "green", T = "red", C = "blue", G = "yellow")
  
  diffList <- list()
  for (i in 2:length(sequences)) {
    diffList[[i - 1]] <- compareSequences(sequences[[1]], sequences[[i]])
    # Map the substitution types to colors
    diffList[[i - 1]]$color <- colorMap[diffList[[i - 1]]$subsType]
    # Default color for other types of substitutions or gaps
    isOther <- is.na(diffList[[i - 1]]$color)
    diffList[[i - 1]]$color[isOther] <- "black"
  }
  
  # Adjust left margin to ensure y-axis labels are not truncated
  oldPar <- par(mar = c(5, 8, 4, 2) + 0.1) # Increase the left margin
  # plot
  plot(NULL, xlim = c(1, genomeLength), ylim = c(0.5, length(sequences)), type = 'n',
       xlab = paste0("Genome Position of ", seqNames[length(seqNames)], ", acting as reference"),
       ylab = "", 
       yaxt = 'n')
  axis(2, at = 1:length(sequences), labels = seqNames, las = 2, cex.axis = 0.8)
  
  # Plot small vertical bars for each difference using mapped colors
  for (i in seq_along(diffList)) {
    yPosStart <- rep(i, nrow(diffList[[i]]))
    yPosEnd <- yPosStart
    for (j in 1:nrow(diffList[[i]])) {
      segments(x0 = diffList[[i]]$position[j], y0 = yPosStart[j] - 0.25, 
               x1 = diffList[[i]]$position[j], y1 = yPosEnd[j] + 0.25, col = diffList[[i]]$color[j])
    }
  }
  
  if (showLegend){
    # Add a semi-transparent legend in the top-left corner
    legend("topleft", inset=c(0.8, 0),
           legend=c("A", "T", "C",  "G", "Other"),
           fill=c( "green", "red", "blue", "yellow", "black"), cex=0.5, bty="n", 
           box.col="gray", bg=adjustcolor("white", alpha.f=0.7))
  }
  
  par(oldPar) # Reset to old graphical parameters after plotting
}



# Example usage:

fastaData <- readFasta(fastaFile = "data/test.fasta", desiredLength = 480)

SNPeek(fastaData, showLegend = F)
