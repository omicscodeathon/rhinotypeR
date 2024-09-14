test_that("pairwiseDistances function works correctly", {
  
  # read a fasta file object for testing
  data_path <- system.file("extdata", "test.fasta", package = "rhinotypeR")
  dna_seqs <- Biostrings::readDNAMultipleAlignment(data_path)
  
  # Determine the number of sequences
  num_sequences <- dim(dna_seqs)[1]
  
  # Test 
  snps_result <- countSNPs(dna_seqs)
  expect_true(is.matrix(snps_result))
  expect_equal(dim(snps_result), c(num_sequences, num_sequences))
  
})
