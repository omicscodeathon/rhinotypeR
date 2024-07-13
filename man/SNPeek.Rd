\name{SNPeek}
\alias{SNPeek}
\title{Visualize single nucleotide polymorphisms}
\description{
Visualizes single nucleotide polymorphisms (SNPs) 
relative to a specified reference sequence. 
To specify the reference, manually move it to the bottom of the alignment. 
Substitutions are color-coded by nucleotide: A = green, T = red, C = blue, G = yellow.
}
\usage{
SNPeek(fastaFile, showLegend = FALSE)
}
\arguments{
  \item{fastaFile}{
  The fasta file used here is the output from the function \code{readFasta}.
  }
  \item{showLegend}{
  Logical indicating whether to show the legend. Default is FALSE.
  }
}
\details{
This function visualizes SNPs in the provided sequence data, 
using a color-coding scheme for different nucleotides.
}
\value{
A plot showing the SNPs relative to a user specified reference sequence.
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

fastaData <- readFasta(fastaFile = test)
SNPeek(fastaData, showLegend = FALSE)
}
\keyword{SNP}
\keyword{visualization}
