
# compare and color sequences based on mutation/ substitutions
compareAndColorSequences <- function(sequences, colorMap, colorFallback = "gray") {
  # Enhanced function to compare sequences and return positions and substitution types
  compareSequences <- function(seqA, seqB) {
    seqAChars <- strsplit(seqA, "")[[1]]
    seqBChars <- strsplit(seqB, "")[[1]]
    differences <- which(seqAChars != seqBChars)
    subsType <- seqBChars[differences]
    data.frame(position = differences, subsType = subsType)
  }
  
  # Compare sequences and assign colors
  diffList <- lapply(2:length(sequences), function(i) {
    diffs <- compareSequences(sequences[[1]], sequences[[i]])
    diffs$color <- colorMap[diffs$subsType]
    diffs$color[is.na(diffs$color)] <- colorFallback
    return(diffs)
  })
  
  return(diffList)
}



# pick respective function based on user's selected model
applyModelFunction <- function(fastaData, model, gapDeletion, functionMap) {
  if (model %in% names(functionMap)) {
    result <- functionMap[[model]](fastaData, gapDeletion = gapDeletion)
  } else {
    stop("Unknown model specified. Choose from 'p-distance', 'JC', 'Kimura2p', or 'Tamura3p' ")
  }
  return(result)
}
