# Global skips for this test file
skip_if_not_installed("MSA2dist")

# Load required data
data(rhinovirusVP4, package = "rhinotypeR")

# Load prototype csv 
ref_path <- system.file("extdata", "prototypes.rda", package = "rhinotypeR")
load(ref_path)  # loads 'prototypes'
prototype_names <- prototypes$Accession

# Shared test fixtures
test_aln_file <- system.file("extdata", "input_aln.fasta", package = "rhinotypeR")
test_aln <- Biostrings::readDNAStringSet(test_aln_file)

# Tests
test_that("assignTypes returns correct data frame structure", {
  skip_if_not(file.exists(test_aln_file), "Test alignment file not found")
  
  result <- assignTypes(test_aln, model = "IUPAC", threshold = 0.105)
  
  expect_s3_class(result, "data.frame")
  expect_named(result, c("query", "assignedType", "distance", "reference"))
  expect_type(result$query, "character")
  expect_type(result$assignedType, "character")
  expect_type(result$distance, "double")
  expect_type(result$reference, "character")
})

test_that("assignTypes errors when prototypes are missing", {
  # Create alignment without prototypes
  no_prototypes <- Biostrings::DNAStringSet(c(
    Query1 = "ATGCATGC",
    Query2 = "ATGGATGC"
  ))
  
  expect_error(
    assignTypes(no_prototypes),
    "prototype sequences"
  )
})

test_that("assignTypes does not include prototypes in output", {
  skip_if_not(file.exists(test_aln_file), "Test alignment file not found")
  
  result <- assignTypes(test_aln, model = "IUPAC", threshold = 0.105)
  
  # Output should not contain any prototype sequences
  expect_false(any(result$query %in% prototype_names))
})

test_that("assignTypes threshold parameter works", {
  skip_if_not(file.exists(test_aln_file), "Test alignment file not found")
  
  # Strict threshold should have more unassigned
  strict <- assignTypes(test_aln, threshold = 0.05)
  lenient <- assignTypes(test_aln, threshold = 0.15)
  
  strict_unassigned <- sum(strict$assignedType == "unassigned")
  lenient_unassigned <- sum(lenient$assignedType == "unassigned")
  
  expect_gte(strict_unassigned, lenient_unassigned)
})

test_that("assignTypes reports distance even when unassigned", {
  skip_if_not(file.exists(test_aln_file), "Test alignment file not found")
  
  # Use very strict threshold
  result <- assignTypes(test_aln, threshold = 0.01)
  
  unassigned_rows <- result[result$assignedType == "unassigned", ]
  
  if (nrow(unassigned_rows) > 0) {
    # Even unassigned should have distance and reference
    expect_false(any(is.na(unassigned_rows$distance)))
    expect_false(any(is.na(unassigned_rows$reference)))
  }
})

test_that("assignTypes model parameter works", {
  skip_if_not(file.exists(test_aln_file), "Test alignment file not found")
  
  expect_no_error(assignTypes(test_aln, model = "IUPAC"))
  expect_no_error(assignTypes(test_aln, model = "raw"))
  expect_no_error(assignTypes(test_aln, model = "JC69"))
})

test_that("assignTypes returns one row per query sequence", {
  skip_if_not(file.exists(test_aln_file), "Test alignment file not found")
  
  result <- assignTypes(test_aln, threshold = 0.105)
  
  # No duplicate query names
  expect_equal(length(unique(result$query)), nrow(result))
})

test_that("assignTypes distance matches reference", {
  skip_if_not(file.exists(test_aln_file), "Test alignment file not found")
  
  result <- assignTypes(test_aln, threshold = 0.105)
  
  # Manually compute distances to verify
  distances <- pairwiseDistances(test_aln, model = "IUPAC")
  
  # Check first query
  first_query <- result$query[1]
  first_ref <- result$reference[1]
  first_dist <- result$distance[1]
  
  if (!is.na(first_ref) && first_query %in% rownames(distances) && 
      first_ref %in% colnames(distances)) {
    expected_dist <- distances[first_query, first_ref]
    expect_equal(first_dist, expected_dist, tolerance = 1e-10)
  }
})
