\name{overallMeanDistance}
\alias{overallMeanDistance}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Calculates the overall mean genetic distance of query sequences using a user-specified evolutionary model
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
overallMeanDistance(RVPrototype, model)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{RVPrototype}{
  %%     ~~Describe \code{pathToRef} here~~
  RVPrototype ids the name of the Rhinovirus prototype. Default is set to "RVAPrototype"
  }
  \item{model}{
  %%     ~~Describe \code{pathToRef} here~~
  It is used to select the model used to assign genotypes to the query sequence. Default is set to p-distance.
  }
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
Martha Luka, Winfred Gatua, Wafaa, Ruth Nanjala
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
  \dontrun{
  # Load the dataset
  test <- system.file("extdata", "input_aln.fasta", package = "rhinotypeR")
  
  # usage
  fastaData <- readFasta(test)
  
  overallMeanDistance(fastaData, model="p-distance")
  overallMeanDistance(fastaData, model="JC")
  overallMeanDistance(fastaData, model="Kimura2p")
  overallMeanDistance(fastaData, model="Tamura3p")
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
