\name{readFasta}
\alias{readFasta}
\title{Reads sequence alignment/fasta files into R for processing}
\description{
Reads a sequence alignment or fasta file into R for processing.
}
\usage{
readFasta(fastaFile)
}
\arguments{
  \item{fastaFile}{
  Path to the fasta file. 
  A FASTA file is a text-based method for encoding nucleotide or peptide sequences, 
  utilizing single-letter codes to represent base pairs or amino acids. 
  Each sequence in FASTA format commences with a header line, 
  succeeded by lines containing the sequence data.
  }
}
\details{
This function reads a FASTA file into R and returns the sequence data for further analysis.
}
\value{
A list containing the sequence data from the FASTA file.
}
\author{
Martha Luka, Ruth Nanjala, Wafaa Rashed, Winfred Gatua, Olaitan Awe 
}
\seealso{
\code{\link{countSNPs}}, \code{\link{SNPeek}}
}
\examples{
# Load the dataset
test <- system.file("extdata", "test.fasta", package = "rhinotypeR")

# Run the command
readFasta(test)
}
\keyword{sequence analysis}
\keyword{data input}
