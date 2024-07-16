\name{plotAA}
\alias{plotAA}
\title{Visualize amino acid substitutions}
\description{
Plots amino acid substitutions with a specified sequence as the reference. 
The input is an amino acid fasta file (translated DNA sequences). 
To specify the reference sequence, move it to the bottom of the alignment. 
Changes are colored by the class of amino acid: 
Red = Positively charged, Blue = Negatively charged, Green = Polar, Yellow = Non-polar.
}
\usage{
plotAA(AAfastaFile, showLegend = FALSE)
}
\arguments{
  \item{AAfastaFile}{
  The file path to the input amino acid sequences in fasta format.
  }
  \item{showLegend}{
  Logical indicating whether to show the legend. Default is FALSE.
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
\code{\link{SNPeek}},
\code{\link{plotAA}}
}
\examples{
# Load the dataset
test <- system.file("extdata", "test.translated.fasta", package = "rhinotypeR")

# Usage
plotAA(test)
}
\keyword{genotype}
\keyword{sequence analysis}
