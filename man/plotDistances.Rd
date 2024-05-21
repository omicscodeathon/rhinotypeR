\name{plotDistances}
\alias{plotDistances}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Visualizes pairwise genetic distances in a heatmap. This function uses the output of pairwiseDistances() as input.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
plotDistances(prototypefile, model = "p-distance")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{prototypefile}{
  %%     ~~Describe \code{pathToRef} here~~
  pathToRef is path to the Rhinovirus prototype. 
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
}
\author{
%%  ~~who you are~~
Martha Luka, Ruth Nanjala, Wafaa Rashed, Winfred Gatua, Olaitan Awe 
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
  
  # Example usage
  fastaD <- readFasta(test)
  distancesToPrototypes <- pairwiseDistances(fastaD, "p-distance")
  plotDistances(distancesToPrototypes)
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
