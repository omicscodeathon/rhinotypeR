\name{plotFrequency}
\alias{plotFrequency}
\title{Plots the frequency of assigned genotypes}
\description{
Plots the frequency of assigned genotypes. 
This function uses the output of \code{assignTypes()} as input.
}
\usage{
plotFrequency(assignedTypesDF, showLegend = FALSE)
}
\arguments{
  \item{assignedTypesDF}{
  A dataframe from the function \code{assignTypes}.
  }
  \item{showLegend}{
  Logical indicating whether to show the legend. Default is FALSE.
  }
}
\details{
This function visualizes the frequency of assigned genotypes based of the newly generated data.
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
