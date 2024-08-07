\name{getPrototypeSeqs}
\alias{getPrototypeSeqs}
\title{Download rhinovirus prototype strains}
\description{
Download rhinovirus prototype strains into a local directory. 
These sequences should be combined with and aligned alongside 
newly generated sequences before being imported into R for genotype assignment.
}
\usage{
getPrototypeSeqs(destinationFolder)
}
\arguments{
  \item{destinationFolder}{
  Provide a path that will act as an output folder for the prototype sequences. 
  The default is set to "output", but you will need to create it 
  prior to running the command in your current directory.
  }
}
\details{
This function downloads rhinovirus prototype strains and saves them in the specified folder.
}
\value{
A character vector of the downloaded file paths.
}
\author{
Martha Luka, Ruth Nanjala, Wafaa Rashed, Winfred Gatua, Olaitan Awe 
}
\seealso{
\code{\link{assignTypes}}
}
\examples{
# Run the function
# getPrototypeSeqs(destinationFolder = "~/Desktop")
}
\keyword{genotype}
\keyword{sequence analysis}
