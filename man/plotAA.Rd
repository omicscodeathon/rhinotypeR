\name{plotAA}
\alias{plotAA}
\title{Visualize amino acid substitutions}
\description{
Plots amino acid substitutions with a specified sequence as the reference. The input is an amino acid fasta file (translated DNA sequences). To specify the reference sequence, move it to the bottom of the alignment. Changes are colored by the class of amino acid: Red = Positively charged, Blue = Negatively charged, Green = Polar, Yellow = Non-polar.
}
\usage{
plotAA(x)
}
\arguments{
  \item{x}{
  The input amino acid sequences.
  }
}
\details{
This function visualizes amino acid substitutions in a given set of sequences with color-coded classes.
}
\value{
A plot object showing the amino acid substitutions.
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

# Usage
fastaData <- readFasta(test)
plotAA(fastaData)
}
\keyword{genotype}
\keyword{sequence analysis}
