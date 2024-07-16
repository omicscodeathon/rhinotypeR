\name{countSNPs}
\alias{countSNPs}
\title{Counts single nucleotide polymorphisms}
\description{
Counts single nucleotide polymorphisms across input sequences and 
produces an output matrix.
}
\usage{
countSNPs(fastaData, gapDeletion = TRUE)
}
\arguments{
  \item{fastaData}{
  The fasta data to be processed.
  }
  \item{gapDeletion}{
  Whether or not to delete positions with gaps. Default is set to TRUE.
  }
}
\details{
This function counts the number of single nucleotide polymorphisms 
(SNPs) across the provided sequences.
}
\value{
A matrix with the SNP counts for each sequence.
}
\author{
Martha Luka, Ruth Nanjala, Wafaa Rashed, Winfred Gatua, Olaitan Awe 
}
\seealso{
\code{\link{readFasta}}
}
\examples{
# Load the dataset
test <- system.file("extdata", "test.fasta", package = "rhinotypeR")

# Run the function
fastaData <- readFasta(test)
countSNPs(fastaData)
}
\keyword{genotype}
\keyword{sequence analysis}
