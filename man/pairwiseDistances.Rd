\name{pairwiseDistances}
\alias{pairwiseDistances}
\title{Estimates pairwise distances}
\description{
Estimates pairwise distances across input sequences using a specified evolutionary model.
}
\usage{
pairwiseDistances(fastaData, model = "p-distance", gapDeletion = TRUE)
}
\arguments{
\item{fastaData}{
  The fasta data to be processed.
  }
 \item{model}{
  The evolutionary model used to calculate distances. Default is set to "p-distance".
  }
  \item{gapDeletion}{
  Whether or not to delete positions with gaps. Default is set to TRUE.
  }
}
\details{
This function calculates the pairwise genetic distances between sequences using the specified evolutionary model.
}
\value{
A matrix of pairwise distances.
}
\author{
Martha Luka, Ruth Nanjala, Wafaa Rashed, Winfred Gatua, Olaitan Awe 
}
\seealso{
\code{\link{overallMeanDistance}}
}
\examples{
# Load the dataset
test <- system.file("extdata", "input_aln.fasta", package = "rhinotypeR")

# Example usage
fastaD <- Biostrings::readDNAStringSet(test)
pairwiseDistances(fastaData = fastaD, model = "p-distance", gapDeletion = TRUE)
}
\keyword{genotype}
\keyword{sequence analysis}
