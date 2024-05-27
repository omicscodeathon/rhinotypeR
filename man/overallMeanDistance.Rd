\name{overallMeanDistance}
\alias{overallMeanDistance}
\title{Estimates the overall mean distance}
\description{
Estimates the overall mean distance of input sequences.
}
\usage{
overallMeanDistance(fastaData, model="p-distance", gapDeletion=TRUE)
}
\arguments{
  \item{fastaFile}{
  The fasta file used here is the output from the function \code{readFasta}.
  }
  \item{model}{
  The evolutionary model used to calculate distances. Default is set to "p-distance".
  }
  \item{gapDeletion}{
  Whether or not to delete positions with gaps. Default is set to TRUE.
  }
}
\details{
This function estimates the overall mean genetic distance of input sequences using the specified evolutionary model.
}
\value{
A numeric value representing the overall mean distance.
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

# Usage
fastaData <- readFasta(test)
overallMeanDistance(fastaData, model="p-distance")
}
\keyword{genotype}
\keyword{sequence analysis}
