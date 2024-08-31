test_that("SNPeek function works correctly", {
  
  # Prepare a small example dataset for testing
  long_sequences <- c("ATGCGCGTATGCGCGT", "ATGCGCGTATGCGCGA", "ATGCGCGTATGCGCGT")
  headers <- c("Seq1", "Seq2", "Seq3")
  fasta_data_large <- list(sequences = long_sequences, headers = headers)
  expect_silent(SNPeek(fasta_data_large, showLegend = FALSE))
  
  
  # Edge case with identical sequences (no differences should be plotted)
  identical_sequences <- c("ATGC", "ATGC", "ATGC")
  fasta_data_identical <- list(sequences = identical_sequences, headers = headers)
  expect_silent(SNPeek(fasta_data_identical, showLegend = FALSE))
  
})
