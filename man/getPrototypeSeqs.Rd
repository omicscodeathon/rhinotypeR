\name{getPrototypeSeqs}
\alias{getPrototypeSeqs}
\title{Download rhinovirus prototype strains}
\description{
Download rhinovirus prototype strains into a local directory. 
These sequences should be combined with and aligned alongside 
newly generated sequences before being imported into R for genotype assignment.
}
\usage{
getPrototypeSeqs(destinationFolder, overwrite = TRUE)
}
\arguments{
  \item{destinationFolder}{
  Provide a path that will act as an output folder for the prototype sequences. 
  The default is set to "output", but you will need to create it 
  prior to running the command in your current directory.
  }
  \item{overwrite}{
  Logical value indicating whether to overwrite existing files in the destination folder. 
  The default is \code{TRUE}.
  }
}
\details{
This function downloads rhinovirus prototype strains and saves them in the specified folder.
If \code{overwrite = TRUE}, existing files with the same name in the destination folder 
will be replaced.
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
if (interactive()) {
    # Specify a destination directory
    dest_dir <- tempdir()  # Create a temporary directory for the example

    # Download prototypes with default overwrite = TRUE
    getPrototypeSeqs(destinationFolder = dest_dir)
    
    # List downloaded files
    list.files(dest_dir)
    }
}
\keyword{genotype}
\keyword{sequence analysis}
