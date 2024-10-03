\name{plotTree}
\alias{plotTree}
\title{Plots a simple phylogenetic tree}
\description{
Plots a simple phylogenetic tree using the genetic distances estimated by 
\code{pairwiseDistances()} and allows for customization using additional arguments.
}
\usage{
plotTree(distance_matrix, ...)
}
\arguments{
  \item{distance_matrix}{
  Distance matrix from the function \code{pairwiseDistances}.
  }
  \item{...}{
  Additional parameters passed to the \code{plot()} function for further customization of the tree plot.
  }
}
\details{
This function visualizes a phylogenetic tree based on the calculated pairwise genetic distances. Users can customize the appearance of the plot by providing additional parameters through the \code{...} argument, which are passed directly to the \code{plot()} function.
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
fastaD <- Biostrings::readDNAStringSet(test)
pdistances <- pairwiseDistances(fastaD, "p-distance")
plotTree(pdistances, col = "blue", lwd = 2)
}
\keyword{genotype}
\keyword{phylogenetics}
