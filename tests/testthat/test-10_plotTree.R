
# Silence plotting devices in tests:
old_dev <- grDevices::dev.cur()
grDevices::pdf(NULL)
on.exit(grDevices::dev.off(), add = TRUE)

testthat::local_edition(3)

# Shared dummy DNAStringSet
simple_dna <- Biostrings::DNAStringSet(c(
  Seq1 = "ATGCATGCATGCATG",
  Seq2 = "ATGCATGCATGCGTG",
  Seq3 = "ATGCATGCATGCTTG",
  Seq4 = "ATGCATGCATGCATG"
))

# Compute distances once
dmat <- pairwiseDistances(simple_dna, model = "IUPAC")

test_that("plotTree returns an hclust object invisibly", {
  hc <- plotTree(dmat)  # will also plot
  expect_s3_class(hc, "hclust")
})

test_that("plotTree accepts matrix or dist objects", {
  hc1 <- plotTree(as.matrix(dmat))
  hc2 <- plotTree(as.dist(dmat))
  expect_s3_class(hc1, "hclust")
  expect_s3_class(hc2, "hclust")
})

test_that("plotTree errors with non-numeric input", {
  bad_matrix <- matrix(letters[1:9], nrow = 3)
  expect_error(plotTree(bad_matrix), "numeric matrix")
})

test_that("plotTree produces symmetric clustering", {
  hc <- plotTree(dmat)
  expect_equal(hc$method, "complete")
  # clustering should have length equal to n - 1 merges
  expect_equal(nrow(hc$merge), length(simple_dna) - 1)
})
