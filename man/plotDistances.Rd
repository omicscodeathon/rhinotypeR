\name{plotDistances}
\alias{plotDistances}
\title{Visualizes pairwise genetic distances}
\description{
Visualizes pairwise genetic distances in a heatmap. This function uses the output of pairwiseDistances() as input.
}
\usage{
plotDistances(distancesMatrix)
}
\arguments{
  \item{distancesMatrix}{
  A matrix of pairwise genetic distances from the function \code{pairwiseDistances}.
  }
}
\details{
This function creates a heatmap to visualize the pairwise genetic distances between sequences.
}
\value{
A heatmap plot object.
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
distancesMatrix <- pairwiseDistances(fastaD, "p-distance")
plotDistances(distancesMatrix)
}
\keyword{genotype}
\keyword{sequence analysis}
