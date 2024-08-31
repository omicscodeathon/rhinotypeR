test_that("pairwiseDistances function works correctly", {
  
  # read a fasta file object for testing
  data_path <- system.file("extdata", "test.fasta", package = "rhinotypeR")
  dna_seqs <- Biostrings::readDNAMultipleAlignment(data_path)
  
  # Determine the number of sequences
  num_sequences <- dim(dna_seqs)[1]
  
  # Test with the default model ("p-distance")
  p_distance_result <- pairwiseDistances(dna_seqs, model = "p-distance")
  expect_true(is.matrix(p_distance_result))
  expect_equal(dim(p_distance_result), c(num_sequences, num_sequences))
  
  # Test with the Jukes-Cantor model ("JC")
  jc_distance_result <- pairwiseDistances(dna_seqs, model = "JC")
  expect_true(is.matrix(jc_distance_result))
  expect_equal(dim(jc_distance_result), c(num_sequences, num_sequences))
  
  # Test with the Kimura 2-parameter model ("Kimura2p")
  k2p_distance_result <- pairwiseDistances(dna_seqs, model = "Kimura2p")
  expect_true(is.matrix(k2p_distance_result))
  expect_equal(dim(k2p_distance_result), c(num_sequences, num_sequences))
  
  # Test with the Tamura 3-parameter model ("Tamura3p")
  t3p_distance_result <- pairwiseDistances(dna_seqs, model = "Tamura3p")
  expect_true(is.matrix(t3p_distance_result))
  expect_equal(dim(t3p_distance_result), c(num_sequences, num_sequences))
  
  # Test with an invalid model
  expect_error(pairwiseDistances(dna_seqs, model = "invalid"), 
               "Unknown model specified. Choose from 'p-distance', 'JC', 'Kimura2p', or 'Tamura3p' ")
})
