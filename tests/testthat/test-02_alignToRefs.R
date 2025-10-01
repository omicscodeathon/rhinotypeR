
# Global skips for this test file
skip_if_not_installed("msa")
skip_if_not_installed("Biostrings") 

# Shared fixture: test fasta path + skip if missing
test_fa <- system.file("extdata", "test.fasta", package = "rhinotypeR")
skip_if_not(file.exists(test_fa), "test.fasta not found in extdata")

# Load sequences once
seqs <- Biostrings::readDNAStringSet(test_fa)

test_that("alignToRefs: basic alignment returns DNAStringSet with equal widths", {
  aln <- alignToRefs(
    seqData   = seqs,
    method    = "ClustalW",
    trimToRef = TRUE,
    refName   = "JN855971.1_A107"
  )
  
  expect_s4_class(aln, "DNAStringSet")
  w <- Biostrings::width(aln)
  expect_true(length(unique(w)) == 1L)
  expect_gt(length(aln), length(seqs))  # prototypes appended
})

test_that("alignToRefs: trimToRef reduces or equals width vs full alignment", {
  aln_full <- alignToRefs(
    seqData   = seqs,
    method    = "ClustalW",
    trimToRef = FALSE
  )
  aln_trim <- alignToRefs(
    seqData   = seqs,
    method    = "ClustalW",
    trimToRef = TRUE,
    refName   = "JN855971.1_A107"
  )
  
  expect_true(all(Biostrings::width(aln_trim) <= Biostrings::width(aln_full)))
})

test_that("alignToRefs: invalid method is rejected", {
  expect_error(
    alignToRefs(seqs, method = "BogusMethod"),
    regexp = "arg should be one of"
  )
})

test_that("alignToRefs: unknown refName errors when trimToRef = TRUE", {
  expect_error(
    alignToRefs(
      seqData   = seqs,
      method    = "ClustalW",
      trimToRef = TRUE,
      refName   = "NOT_IN_ALIGNMENT_XXX"
    ),
    regexp = "refName not found|Reference name not found"
  )
})

test_that("alignToRefs: result is suitable for assignTypes (contains prototypes)", {
  res <- assignTypes(
    alignToRefs(
      seqData   = seqs,
      method    = "ClustalW",
      trimToRef = TRUE,
      refName   = "JN855971.1_A107"
    ),
    model = "IUPAC",
    threshold = 0.105
  )
  
  expect_s3_class(res, "data.frame")
  expect_true(all(c("query","assignedType","distance","reference") %in% names(res)))
})
