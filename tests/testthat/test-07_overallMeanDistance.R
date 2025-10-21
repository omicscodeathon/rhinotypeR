
# Global skips for this test file
skip_if_not_installed("MSA2dist")

# Small deterministic 15-nt alignment (no gaps)
simple_dna <- Biostrings::DNAStringSet(c(
  Seq1 = "ATGCATGCATGCATG",
  Seq2 = "ATGGATGCATGCATG",  # 1 mismatch vs Seq1
  Seq3 = "ATGCTTGCATGCATG",  # a couple of mismatches vs Seq1
  Seq4 = "ATGCATGCATGCATG"   # identical to Seq1
))

# Same length (15) but with gaps in some sequences
# (global deletion should remove that column)
simple_dna_gaps <- Biostrings::DNAStringSet(c(
  Seq1 = "ATGCA-TGCATGCAT",  
  Seq2 = "ATGCAATGCA-GCAT",
  Seq3 = "ATGCAATGCATGCAT",
  Seq4 = "AT-CAATGCATGCAT"
))

test_that("overallMeanDistance returns a single finite numeric in [0,1]", {
  m <- overallMeanDistance(simple_dna, model = "IUPAC", deleteGapsGlobally = FALSE)
  expect_type(m, "double")
  expect_length(m, 1L)
  expect_true(is.finite(m))
  expect_gte(m, 0)
  expect_lte(m, 1)
})

test_that("overallMeanDistance agrees with mean(lower triangle) of pairwiseDistances", {
  m_fun <- overallMeanDistance(simple_dna, model = "IUPAC", deleteGapsGlobally = FALSE)
  
  D <- pairwiseDistances(simple_dna, model = "IUPAC", deleteGapsGlobally = FALSE)
  Dm <- as.matrix(D)
  m_manual <- mean(Dm[lower.tri(Dm)])
  
  expect_equal(m_fun, m_manual, tolerance = 1e-12)
})

test_that("deleteGapsGlobally does not increase mean distance", {
  m_keep <- overallMeanDistance(simple_dna_gaps, model = "IUPAC", deleteGapsGlobally = FALSE)
  m_drop <- overallMeanDistance(simple_dna_gaps, model = "IUPAC", deleteGapsGlobally = TRUE)
  
  # Removing a gap column across all sequences should not increase mean distance
  expect_lte(m_drop, m_keep)
})

test_that("overallMeanDistance runs with common models", {
  expect_no_error(overallMeanDistance(simple_dna, model = "IUPAC"))
  expect_no_error(overallMeanDistance(simple_dna, model = "raw"))
  expect_no_error(overallMeanDistance(simple_dna, model = "JC69"))
})

test_that("Two-sequence case equals their pairwise distance (raw model)", {
  x <- Biostrings::DNAStringSet(c(
    S1 = "ATGCATGCATGCATG",
    S2 = "ATGCATGCATGCATT"  
  ))
  
  m <- overallMeanDistance(x, model = "raw", deleteGapsGlobally = FALSE)
  
  D <- pairwiseDistances(x, model = "raw", deleteGapsGlobally = FALSE)
  expected <- as.matrix(D)["S1", "S2"]
  
  expect_equal(m, expected, tolerance = 1e-12)
})

test_that("Identical sequences yield mean distance 0", {
  x <- Biostrings::DNAStringSet(c(
    A = "ACGTACGTACGTACG",
    B = "ACGTACGTACGTACG",
    C = "ACGTACGTACGTACG"
  ))
  m <- overallMeanDistance(x, model = "IUPAC")
  expect_equal(m, 0)
})

test_that("Non-DNAStringSet input errors (type safety)", {
  expect_error(overallMeanDistance(list(A = "ACGT")), regexp = "")
})
