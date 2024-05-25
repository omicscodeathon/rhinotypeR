\name{pairwiseDistances}
\alias{pairwiseDistances}
\title{Estimates pairwise distances}
\description{
Estimates pairwise distances across input sequences using a specified evolutionary model.
}
\usage{
pairwiseDistances(prototypefile, model = "p-distance")
}
\arguments{
\item{prototypefile}{
  The path to the Rhinovirus prototype sequences.
  }
 \item{model}{
  The evolutionary model used to calculate distances. Default is set to "p-distance".
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
fastaD <- readFasta(test)
pairwiseDistances(fastaD, "p-distance")
}
\keyword{genotype}
\keyword{sequence analysis}
