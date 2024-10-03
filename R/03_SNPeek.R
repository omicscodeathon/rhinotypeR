SNPeek <- function(fastaData, showLegend = FALSE) {
  
  # Preprocess fasta data
  fastaData <- preProcessFastaStringSet(fastaData)
  
  sequences <- fastaData$sequences
  seqNames <- fastaData$headers
  genomeLength <- max(vapply(sequences, nchar, integer(1)))
  
  # Define a color map for nucleotides
  colorMapNT <- c(A = "green", T = "red", C = "blue", G = "yellow")
  
  # Use helper function to compare sequences and assign colors
  diffList <- compareAndColorSequences(sequences, colorMapNT, colorFallback = "black")
  
  # Adjust left margin to ensure y-axis labels are not truncated
  oldPar <- par(mar = c(5, 8, 4, 2) + 0.1)
  
  plot(NULL, xlim = c(1, genomeLength), ylim = c(0.5, length(sequences)), type = 'n',
       xlab = paste0("Genome Position of ", seqNames[length(seqNames)], ", acting as reference"),
       ylab = "", yaxt = 'n')
  
  axis(2, at = seq_along(sequences), labels = seqNames, las = 2, cex.axis = 0.8)
  
  # Plot small vertical bars for each difference using mapply
  mapply(function(diffs, i) {
    segments(x0 = diffs$position, y0 = i - 0.25, 
             x1 = diffs$position, y1 = i + 0.25, col = diffs$color)
  }, diffList, seq_along(diffList))
  
  if (showLegend) {
    legend("topleft", inset = c(0.8, 0),
           legend = c("A", "T", "C", "G", "Other"),
           fill = c("green", "red", "blue", "yellow", "black"), cex = 0.5, bty = "n", 
           box.col = "gray", bg = adjustcolor("white", alpha.f = 0.7))
  }
  
  par(oldPar) # Reset to old graphical parameters after plotting
}
