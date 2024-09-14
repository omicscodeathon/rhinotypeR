\name{assignTypes}
\alias{assignTypes}
\title{Assigns genotypes to input sequences}
\description{
Rapidly assigns genotypes to input sequences. 
The input fasta file should include the prototype strains, 
which can be downloaded using getPrototypeSeqs().
}
\usage{
assignTypes(fastaData, model = "p-distance", gapDeletion = TRUE, threshold = 0.105)
}
\arguments{
 \item{fastaData}{
  The fasta data to be processed.
  }
  \item{model}{
  The evolutionary model to be used. Default is set to "p-distance".
  }
  \item{gapDeletion}{
  Logical. If TRUE, gaps are deleted. Default is TRUE.
  }
  \item{threshold}{
  The distance threshold for genotype assignment. Default is 0.105.
  }
}
\details{
This function compares input sequences against prototype strains using a 
specified evolutionary model and assigns genotypes based on predefined distance thresholds.
}
\value{
A list with the assigned genotypes and their distances to the reference sequences.
}
\author{
Martha Luka, Ruth Nanjala, Wafaa Rashed, Winfred Gatua, Olaitan Awe 
}
\note{
Ensure the input fasta file includes the necessary prototype strains.
}
\seealso{
\code{\link{getPrototypeSeqs}}, \code{\link{pairwiseDistances}}
}
\examples{
# Load the dataset
test <- system.file("extdata", "input_aln.fasta", package = "rhinotypeR")

# Run command
fastaD <- Biostrings::readDNAStringSet(test)
assignTypes(fastaD, model = "p-distance", gapDeletion = TRUE, threshold = 0.105)
}
\keyword{genotype}
\keyword{sequence analysis}
