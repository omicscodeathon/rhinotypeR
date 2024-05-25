\name{plotTree}
\alias{plotTree}
\title{Plots a simple phylogenetic tree}
\description{
Plots a simple phylogenetic tree using the genetic distances estimated by \code{pairwiseDistances()}.
}
\usage{
plotTree(distance_matrix)
}
\arguments{
  \item{distance_matrix}{
  Distance matrix from the function \code{pairwiseDistances}.
  }
}
\details{
This function visualizes a phylogenetic tree based on the calculated pairwise genetic distances.
}
\value{
A plot object representing the phylogenetic tree.
}
\author{
Martha Luka, Ruth Nanjala, Wafaa Rashed, Winfred Gatua, Olaitan Awe 
}
\seealso{
\code{\link{pairwiseDistances}}
}
\examples{
# Load the dataset
test <- system.file("extdata", "input_aln.fasta", package = "rhinotypeR")

# Example usage
fastaD <- readFasta(test)
pdistances <- pairwiseDistances(fastaD, "p-distance")
plotTree(pdistances)
}
\keyword{genotype}
\keyword{phylogenetics}
