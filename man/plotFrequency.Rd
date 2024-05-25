\name{plotFrequency}
\alias{plotFrequency}
\title{Plots the frequency of assigned genotypes}
\description{
Plots the frequency of assigned genotypes. This function uses the output of \code{assignTypes()} as input.
}
\usage{
plotFrequency(target_data, ref_data, model = "Tamura3p")
}
\arguments{
  \item{target_data}{
  Path to the target fasta file.
  }
  \item{ref_data}{
  Path to the reference fasta file.
  }
  \item{model}{
  The evolutionary model used for genotype assignment. Default is "Tamura3p".
  }
}
\details{
This function visualizes the frequency of assigned genotypes based on the provided target and reference data.
}
\value{
A plot showing the frequency of assigned genotypes.
}
\author{
Martha Luka, Ruth Nanjala, Wafaa Rashed, Winfred Gatua, Olaitan Awe 
}
\seealso{
\code{\link{assignTypes}}
}
\examples{
# Load the dataset
test <- system.file("extdata", "input_aln.fasta", package = "rhinotypeR")

# Run 
queryFastaData <- readFasta(test)
df <- assignTypes(queryFastaData, "p-distance")

plotFrequency(df)
}
\keyword{genotype}
\keyword{visualization}
