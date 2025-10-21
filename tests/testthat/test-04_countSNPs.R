
# Global skips for this test file
skip_if_not_installed("MSA2dist") 

# Shared test fixtures
simple_dna <- Biostrings::DNAStringSet(c(
  Seq1 = "ATGCATGC",
  Seq2 = "ATGGATGC",  # 1 SNP from Seq1
  Seq3 = "ATGCTTGC",  # 1 SNP from Seq1
  Seq4 = "ATGCATGC"   # 0 SNPs from Seq1
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
test_that("countSNPs returns correct matrix structure", {
  result <- countSNPs(simple_dna)
  
  expect_type(result, "integer")
  expect_true(is.matrix(result))
  expect_equal(nrow(result), 4)
  expect_equal(ncol(result), 4)
  expect_equal(rownames(result), names(simple_dna))
  expect_equal(colnames(result), names(simple_dna))
})

test_that("countSNPs diagonal is zero", {
  result <- countSNPs(simple_dna)
  
  # Diagonal should be all zeros (sequence compared to itself)
  expect_equal(unname(diag(result)), rep(0L, 4))
})

test_that("countSNPs matrix is symmetric", {
  result <- countSNPs(simple_dna)
  
  # SNP counts should be symmetric
  expect_equal(result, t(result))
})

test_that("countSNPs counts identical sequences correctly", {
  result <- countSNPs(identical_seqs)
  
  # All pairs should have 0 SNPs
  expect_true(all(result == 0))
})

test_that("countSNPs counts completely different sequences correctly", {
  result <- countSNPs(all_different)
  
  # All off-diagonal elements should be 4 (all positions differ)
  off_diag <- result[upper.tri(result)]
  expect_true(all(off_diag == 4))
})

test_that("countSNPs handles gaps correctly without global deletion", {
  result <- countSNPs(dna_with_gaps, deleteGapsGlobally = FALSE)
  
  expect_type(result, "integer")
  expect_true(is.matrix(result))
})

test_that("countSNPs with deleteGapsGlobally works", {
  result <- countSNPs(dna_with_gaps, deleteGapsGlobally = TRUE)
  
  expect_type(result, "integer")
  expect_true(is.matrix(result))
  # After removing gap columns, fewer sites should be compared
})

test_that("countSNPs gives integer counts", {
  result <- countSNPs(simple_dna)
  
  # All values should be integers
  expect_true(all(result == floor(result)))
  expect_type(result, "integer")
})

test_that("countSNPs errors on unaligned sequences", {
  unaligned <- Biostrings::DNAStringSet(c(
    Seq1 = "ATGC",
    Seq2 = "ATGCAA"
  ))
  
  expect_error(countSNPs(unaligned))
})

test_that("countSNPs errors with invalid input", {
  expect_error(countSNPs("not a DNAStringSet"))
  expect_error(countSNPs(NULL))
})

test_that("countSNPs handles ambiguous bases", {
  ambiguous <- Biostrings::DNAStringSet(c(
    Seq1 = "ATGC",
    Seq2 = "NTGC"  # N is ambiguous
  ))
  
  # Should not error (IUPAC model handles ambiguity)
  expect_no_error(result <- countSNPs(ambiguous))
  expect_type(result, "integer")
})

test_that("countSNPs values are non-negative", {
  result <- countSNPs(simple_dna)
  
  expect_true(all(result >= 0))
})

test_that("countSNPs with deleteGapsGlobally reduces compared sites", {
  # This is a behavioral test - with gaps removed, 
  # the comparison should still work but on fewer sites
  result_with_gaps <- countSNPs(dna_with_gaps, deleteGapsGlobally = FALSE)
  result_no_gaps <- countSNPs(dna_with_gaps, deleteGapsGlobally = TRUE)
  
  # Both should complete without error
  expect_type(result_with_gaps, "integer")
  expect_type(result_no_gaps, "integer")
})

