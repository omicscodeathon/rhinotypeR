test_that("assignTypes function works correctly", {
  
  # Read a FASTA file object for testing
  comprehensiveData_path <- system.file("extdata", "input_aln.fasta", package = "rhinotypeR")
  comprehensiveData <- Biostrings::readDNAMultipleAlignment(comprehensiveData_path)
  
  # Test when prototypes are present in the fasta data
  result <- assignTypes(comprehensiveData, model = "p-distance", gapDeletion = TRUE, threshold = 0.105)
  expect_true(is.data.frame(result))
  expect_true(all(c("query", "assignedType", "distance", "reference") %in% colnames(result)))

  # Test with an altered FASTA file where prototypes are missing
  incomprehensiveData_path <- system.file("extdata", "test.fasta", package = "rhinotypeR")
  incomprehensiveData <- Biostrings::readDNAMultipleAlignment(incomprehensiveData_path)
  
  
  expect_error(assignTypes(incomprehensiveData, model = "p-distance", gapDeletion = TRUE, threshold = 0.105),
               "To classify rhinovirus sequences, please ensure your input sequences contain the prototypes.")
  
  # Check that the function assigns types correctly when valid data and threshold are provided
  correct_result <- assignTypes(comprehensiveData, model = "p-distance", gapDeletion = TRUE, threshold = 0.05)
  correct_result <- na.omit(correct_result) # drop NAs which are values above the threshold
  expect_true(all(correct_result$distance <= 0.05))

})

