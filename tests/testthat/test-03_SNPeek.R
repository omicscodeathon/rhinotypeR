
# Silence plotting devices in tests:
testthat::local_null_device()
testthat::local_edition(3)

# Shared test fixtures
simple_aln <- Biostrings::DNAStringSet(c(
  Seq1 = "ATGCATGC",
  Seq2 = "ATGGATGC",
  Seq3 = "ATGCTTGC",
  Ref  = "ATGCATGC"
))

short_aln <- Biostrings::DNAStringSet(c(
  Seq1 = "ATGC",
  Seq2 = "ATGG",
  Seq3 = "ATGC"
))

long_aln <- Biostrings::DNAStringSet(c(
  Seq1 = "ATGCATGCATGC",
  Seq2 = "ATGGATGCATGC",
  Ref  = "ATGCATGCATGC"
))

aln_with_gaps <- Biostrings::DNAStringSet(c(
  Seq1 = "ATG-ATGC",
  Seq2 = "ATGGATGC",
  Ref  = "ATGCATGC"
))

# Tests
test_that("SNPeek works with basic aligned sequences", {
  expect_silent(result <- SNPeek(simple_aln))
  expect_s3_class(result, "SNPeekCache")
  expect_equal(result$genome_len, 8)
})

test_that("SNPeek handles reference selection correctly", {
  result1 <- SNPeek(short_aln)
  expect_equal(result1$ref_label, "Seq3")
  
  result2 <- SNPeek(short_aln, ref_name = "Seq1")
  expect_equal(result2$ref_label, "Seq1")
})

test_that("SNPeek errors on unaligned sequences", {
  unaligned <- Biostrings::DNAStringSet(c(Seq1 = "ATGC", Seq2 = "ATGCAA"))
  expect_error(SNPeek(unaligned), "same width")
})

test_that("SNPeek errors with invalid input type", {
  expect_error(SNPeek("not valid"), "DNAStringSet")
})

test_that("SNPeek xlim and window parameters work", {
  expect_silent(SNPeek(long_aln, xlim = c(3, 8)))
  expect_silent(SNPeek(long_aln, center = 6, window = 4))
})

test_that("SNPeek highlighting works", {
  expect_silent(SNPeek(short_aln, highlight_seqs = c("Seq1", "Seq2")))
  expect_warning(SNPeek(short_aln, highlight_seqs = "NonExistent"), "not found")
})

test_that("SNPeek show_only_highlighted filters correctly", {
  result <- SNPeek(simple_aln, 
                   highlight_seqs = c("Seq1", "Seq2"),
                   show_only_highlighted = TRUE)
  expect_s3_class(result, "SNPeekCache")
})

test_that("SNPeek custom color maps work", {
  custom_colors <- c(A = "pink", T = "purple", C = "orange", G = "cyan")
  result <- SNPeek(short_aln, colorMapNT = custom_colors)
  expect_equal(result$colorMap, custom_colors)
})

test_that("SNPeek handles sequences with gaps", {
  expect_silent(result <- SNPeek(aln_with_gaps))
  expect_s3_class(result, "SNPeekCache")
})

test_that("SNPeek legend parameter works", {
  expect_silent(SNPeek(short_aln, showLegend = TRUE))
  expect_silent(SNPeek(short_aln, showLegend = FALSE))
})

test_that("SNPeek cache structure is correct", {
  cache <- SNPeek(short_aln)
  
  expect_named(cache, c("genome_len", "seq_names", "ref_label", "ref_row",
                        "diffs_list", "col_levels", "colorMap", "colorFallback"))
  expect_type(cache$genome_len, "integer")
  expect_type(cache$seq_names, "character")
  expect_type(cache$diffs_list, "list")
})

test_that("SNPeek works with package example data", {
  skip_if_not_installed("rhinotypeR")
  data(rhinovirusVP4, package = "rhinotypeR")
  
  expect_silent(result <- SNPeek(rhinovirusVP4))
  expect_s3_class(result, "SNPeekCache")
  expect_true(result$genome_len > 0)
})