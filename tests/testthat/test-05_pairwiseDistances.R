
# Global skips for this test file
skip_if_not_installed("MSA2dist")

# Shared test fixtures
simple_dna <- Biostrings::DNAStringSet(c(
  Seq1 = "ATGCATGC",
  Seq2 = "ATGGATGC",
  Seq3 = "ATGCTTGC",
  Seq4 = "ATGCATGC"
))

dna_with_gaps <- Biostrings::DNAStringSet(c(
  Seq1 = "ATGC-TGC",
  Seq2 = "ATGGATGC",
  Seq3 = "ATGC-TGC"
))

identical_seqs <- Biostrings::DNAStringSet(c(
  Seq1 = "ATGC",
  Seq2 = "ATGC",
  Seq3 = "ATGC"
))

all_different <- Biostrings::DNAStringSet(c(
  Seq1 = "AAAA",
  Seq2 = "TTTT",
  Seq3 = "CCCC"
))

# Tests
test_that("pairwiseDistances returns correct matrix structure", {
  result <- pairwiseDistances(simple_dna)
  
  expect_type(result, "double")
  expect_true(is.matrix(result))
  expect_equal(nrow(result), 4)
  expect_equal(ncol(result), 4)
  expect_equal(rownames(result), names(simple_dna))
  expect_equal(colnames(result), names(simple_dna))
})

test_that("pairwiseDistances diagonal is zero", {
  result <- pairwiseDistances(simple_dna)
  
  # Diagonal should be all zeros (sequence compared to itself)
  expect_equal(diag(result), rep(0, 4))
})

test_that("pairwiseDistances matrix is symmetric", {
  result <- pairwiseDistances(simple_dna)
  
  # Distance matrix should be symmetric
  expect_equal(result, t(result))
})

test_that("pairwiseDistances handles identical sequences", {
  result <- pairwiseDistances(identical_seqs)
  
  # All pairs should have 0 distance
  expect_true(all(result == 0))
})

test_that("pairwiseDistances handles completely different sequences", {
  result <- pairwiseDistances(all_different)
  
  # All off-diagonal elements should be > 0
  off_diag <- result[upper.tri(result)]
  expect_true(all(off_diag > 0))
  expect_true(all(off_diag <= 1))  # Distances should be proportions (0-1)
})

test_that("pairwiseDistances values are in valid range", {
  result <- pairwiseDistances(simple_dna)
  
  # All distances should be between 0 and 1
  expect_true(all(result >= 0))
  expect_true(all(result <= 1))
})

test_that("pairwiseDistances model parameter works", {
  # Should accept different models
  expect_silent(pairwiseDistances(simple_dna, model = "IUPAC"))
  expect_silent(pairwiseDistances(simple_dna, model = "raw"))
  expect_silent(pairwiseDistances(simple_dna, model = "JC69"))
})

test_that("pairwiseDistances deleteGapsGlobally works", {
  result_with_gaps <- pairwiseDistances(dna_with_gaps, deleteGapsGlobally = FALSE)
  result_no_gaps <- pairwiseDistances(dna_with_gaps, deleteGapsGlobally = TRUE)
  
  expect_type(result_with_gaps, "double")
  expect_type(result_no_gaps, "double")
  
  # Both should be valid distance matrices
  expect_true(is.matrix(result_with_gaps))
  expect_true(is.matrix(result_no_gaps))
})

test_that("pairwiseDistances errors on unaligned sequences", {
  unaligned <- Biostrings::DNAStringSet(c(
    Seq1 = "ATGC",
    Seq2 = "ATGCAA"
  ))
  
  expect_error(pairwiseDistances(unaligned))
})

test_that("pairwiseDistances errors with invalid input", {
  expect_error(pairwiseDistances("not a DNAStringSet"))
  expect_error(pairwiseDistances(NULL))
})

test_that("pairwiseDistances handles ambiguous bases", {
  ambiguous <- Biostrings::DNAStringSet(c(
    Seq1 = "ATGC",
    Seq2 = "NTGC"
  ))
  
  # IUPAC model should handle ambiguity
  expect_silent(result <- pairwiseDistances(ambiguous, model = "IUPAC"))
  expect_type(result, "double")
})
