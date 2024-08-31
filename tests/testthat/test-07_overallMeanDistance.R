test_that("overallMeanDistance function works correctly", {
  
  # read a fasta file object for testing
  data_path <- system.file("extdata", "test.fasta", package = "rhinotypeR")
  dna_seqs <- Biostrings::readDNAMultipleAlignment(data_path)
  
  # Test with the default model ("p-distance")
  p_distance_result <- overallMeanDistance(dna_seqs, model = "p-distance")
  expect_true(is.numeric(p_distance_result))
  expect_equal(round(p_distance_result, 2), 0.22) # values agreed by other software (MegaX and ape R package)
  
  # Test with the Jukes-Cantor model ("JC")
  jc_distance_result <- overallMeanDistance(dna_seqs, model = "JC")
  expect_true(is.numeric(jc_distance_result))
  expect_equal(round(jc_distance_result, 2), 0.26) # values agreed by other software (MegaX and ape R package)
  
  # Test with the Kimura 2-parameter model ("Kimura2p")
  k2p_distance_result <- overallMeanDistance(dna_seqs, model = "Kimura2p")
  expect_true(is.numeric(k2p_distance_result))
  expect_false(is.na(k2p_distance_result)) # Ensure that a valid result is returned
  
  # Test with the Tamura 3-parameter model ("Tamura3p")
  t3p_distance_result <- overallMeanDistance(dna_seqs, model = "Tamura3p")
  expect_true(is.numeric(t3p_distance_result))
  expect_false(is.na(t3p_distance_result)) # Ensure that a valid result is returned
  
  # Test with an invalid model
  expect_error(overallMeanDistance(dna_seqs, model = "invalid"), 
               "Unknown model specified. Choose from 'p-distance', 'JC', 'Kimura2p', or 'Tamura3p' ")
})