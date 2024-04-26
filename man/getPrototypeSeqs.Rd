\name{getPrototypeSeqs}
\alias{getPrototypeSeqs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Downloads RV prototypes into a user-specified local directory

}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
getPrototypeSeqs(destinationFolder = "./output")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{destinationFolder}{
%%     ~~Describe \code{x} here~~
  Provide a path that will act as an output folder for the Prototype sequences.
  The default is set to "output", but you will need to create it prior to running the command in your current directory.
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
Martha Luka, Winfred Gatua, Wafaa, Ruth Nanjala
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
  getPrototypeSeqs(destinationFolder = "./output")
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
