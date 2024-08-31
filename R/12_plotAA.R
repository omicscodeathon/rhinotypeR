plotAA <- function(AAfastaFile, showLegend = FALSE) {
  
  # Preprocess fasta data
  AAfastaFile <- preProcessFastaStringSet(AAfastaFile)
  
  sequences <- AAfastaFile$sequences
  seqNames <- AAfastaFile$headers
  proteinLength <- max(vapply(sequences, nchar, integer(1)))
  
  # Define a color map for amino acids
  colorMap <- c(
    R = "red", H = "red", K = "red",      # Positively charged amino acid 
    D = "blue", E = "blue",               # Negatively charged amino acid 
    S = "green", T = "green", N = "green", Q = "green", # Polar amino acid 
    A = "yellow", V = "yellow", I = "yellow", L = "yellow", M = "yellow", F = "yellow", 
    W = "yellow", P = "yellow", G = "yellow", Y = "yellow", C = "yellow"  # Nonpolar amino acid 
  )
  
  # Use helper function to compare sequences and assign colors
  diffList <- compareAndColorSequences(sequences, colorMap, colorFallback = "gray")
  
  # Adjust left margin to ensure y-axis labels are not truncated
  oldPar <- par(mar = c(5, 8, 4, 2) + 0.1, xpd = TRUE)
  
  plot(NULL, xlim = c(1, proteinLength), ylim = c(0.5, length(sequences)), type = 'n',
       xlab = paste0("Protein Position of ", seqNames[length(seqNames)], ", acting as reference"),
       ylab = "", yaxt = 'n')
  
  axis(2, at = seq_along(sequences), labels = seqNames, las = 2, cex.axis = 0.8)
  
  # Plot small vertical bars for each difference using mapply
  mapply(function(diffs, i) {
    segments(x0 = diffs$position, y0 = i - 0.25, 
             x1 = diffs$position, y1 = i + 0.25, col = diffs$color)
  }, diffList, seq_along(diffList))
  
  if (showLegend) { 
    legend("topleft", inset = c(0.78, 0),
           legend = c("+ve charged", "-ve charged", "Polar", "Non-polar", "Other"),
           fill = c("red", "blue", "green", "yellow", "gray"), cex = 0.45, bty = "n", 
           box.col = "gray", bg = adjustcolor("white", alpha.f = 0.7))
  }
  
  par(oldPar) # Reset to old graphical parameters after plotting
}
