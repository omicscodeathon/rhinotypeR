source("./scripts/02_readFasta.R")

plotAA <- function(fastaData, showLegend = FALSE) {
  sequences <- fastaData$sequences
  seqNames <- fastaData$headers
  proteinLength <- max(sapply(sequences, nchar))
  
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
  legend("topleft", inset=c(0.78, 0),
         legend=c("+ve charged", "-ve charged", "Polar", "Non-polar", "Other"),
         fill=c("red", "blue", "green", "yellow", "gray"), cex=0.45, bty="n", 
         box.col="gray", bg=adjustcolor("white", alpha.f=0.7))
  
  par(oldPar) # Reset to old graphical parameters after plotting
  }
}
