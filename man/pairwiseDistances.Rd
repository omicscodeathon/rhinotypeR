\name{pairwiseDistances}
\alias{pairwiseDistances}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Calculate pairwise distance among input sequences
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
pairwiseDistances(prototypefile, model = "p-distance")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{RVPrototype}{
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
##---- Should be DIRECTLY executable !! ----

## The function is currently defined as
pairwiseDistances("inputSequencesPath", "p-distance")

## An example on the test dataset will be:
pairwiseDistances(RVBPrototype, "p-distance")
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