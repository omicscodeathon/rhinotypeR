# Silence plotting devices in tests:
testthat::local_null_device()
testthat::local_edition(3)

# Shared test fixtures
simple_aa <- Biostrings::AAStringSet(c(
  Seq1 = "MTEYKLVVVGYKL",
  Seq2 = "MTEYKLVILVVVG",
  Seq3 = "MTEYKLVVV-LVV",
  Ref  = "MTEYKLVVVGYKL"
))

short_aa <- Biostrings::AAStringSet(c(
  Seq1 = "ARNDQ",
  Seq2 = "ARNDE",
  Seq3 = "ARNDQ"
))

long_aa <- Biostrings::AAStringSet(c(
  Seq1 = "ARNDCQEGHILKMFPSTWYV",
  Seq2 = "ARNDCQEGHILKMFPSTWYA",
  Ref  = "ARNDCQEGHILKMFPSTWYV"
))

aa_with_gaps <- Biostrings::AAStringSet(c(
  Seq1 = "ARN-CQEG",
  Seq2 = "ARNDCQEG",
  Ref  = "ARNDCQEG"
))

# Tests
test_that("plotAA works with basic aligned sequences", {
  expect_silent(result <- plotAA(simple_aa))
  expect_s3_class(result, "SNPeekCache")
  expect_equal(result$genome_len, 13)
})

test_that("plotAA handles reference selection correctly", {
  result1 <- plotAA(short_aa)
  expect_equal(result1$ref_label, "Seq3")
  
  result2 <- plotAA(short_aa, ref_name = "Seq1")
  expect_equal(result2$ref_label, "Seq1")
})

test_that("plotAA errors on unaligned sequences", {
  unaligned <- Biostrings::AAStringSet(c(Seq1 = "ARND", Seq2 = "ARNDCQ"))
  expect_error(plotAA(unaligned), "same width")
})

test_that("plotAA errors with invalid input type", {
  expect_error(plotAA("not valid"), "AAStringSet")
  
  # Should also error with DNAStringSet
  dna <- Biostrings::DNAStringSet(c(Seq1 = "ATGC", Seq2 = "ATGG"))
  expect_error(plotAA(dna), "AAStringSet")
})

test_that("plotAA xlim and window parameters work", {
  expect_silent(plotAA(long_aa, xlim = c(5, 15)))
  expect_silent(plotAA(long_aa, center = 10, window = 8))
})

test_that("plotAA highlighting works", {
  expect_silent(plotAA(short_aa, highlight_seqs = c("Seq1", "Seq2")))
  expect_warning(plotAA(short_aa, highlight_seqs = "NonExistent"), "not found")
})

test_that("plotAA show_only_highlighted filters correctly", {
  result <- plotAA(simple_aa, 
                   highlight_seqs = c("Seq1", "Seq2"),
                   show_only_highlighted = TRUE)
  expect_s3_class(result, "SNPeekCache")
})

test_that("plotAA custom color maps work", {
  custom_colors <- c(A = "pink", R = "purple", N = "orange", D = "cyan", Q = "brown")
  result <- plotAA(short_aa, colorMapAA = custom_colors)
  expect_true(all(names(custom_colors) %in% names(result$colorMap)))
})

test_that("plotAA handles sequences with gaps", {
  expect_silent(result <- plotAA(aa_with_gaps))
  expect_s3_class(result, "SNPeekCache")
})

test_that("plotAA legend parameter works", {
  expect_silent(plotAA(short_aa, showLegend = TRUE))
  expect_silent(plotAA(short_aa, showLegend = FALSE))
})

test_that("plotAA default color map covers all standard amino acids", {
  result <- plotAA(long_aa)
  
  # Check that all 20 standard amino acids are in the default map
  standard_aa <- c("A", "R", "N", "D", "C", "Q", "E", "G", "H", "I",
                   "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V")
  expect_true(all(standard_aa %in% names(result$colorMap)))
})

test_that("plotAA cache structure is correct", {
  cache <- plotAA(short_aa)
  
  expect_named(cache, c("genome_len", "seq_names", "ref_label", "ref_row",
                        "diffs_list", "col_levels", "colorMap", "colorFallback"))
  expect_type(cache$genome_len, "integer")
  expect_type(cache$seq_names, "character")
  expect_type(cache$diffs_list, "list")
})

test_that("plotAA handles unknown amino acids", {
  aa_with_x <- Biostrings::AAStringSet(c(
    Seq1 = "ARNDX",
    Ref  = "ARNDC"
  ))
  
  expect_silent(result <- plotAA(aa_with_x))
  expect_s3_class(result, "SNPeekCache")
})