test_that("plotTree function works correctly", {
  
  # read a fasta file object for testing
  data_path <- system.file("extdata", "input_aln.fasta", package = "rhinotypeR")
  dna_seqs <- Biostrings::readDNAMultipleAlignment(data_path)
  
  # Generate a test matrix using the pairwiseDistances function
  distance_matrix <- pairwiseDistances(dna_seqs, model = "p-distance")
  
  # Test that distance matrix is correctly converted to a "dist" object
  dist_object <- as.dist(distance_matrix)
  expect_true(inherits(dist_object, "dist"))
  
  # Test that hierarchical clustering works without error
  hc <- hclust(dist_object, method = "complete")
  expect_true(inherits(hc, "hclust"))
  
  # Test that the plot function works without errors (no visual verification)
  expect_error(plotTree(distance_matrix), NA)  # NA means no error is expected
  
  # Optionally, you can check properties of the hc object, such as the order of labels
  expect_equal(hc$labels, rownames(distance_matrix))
})
