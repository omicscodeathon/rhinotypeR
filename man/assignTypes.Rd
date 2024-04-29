\name{assignTypes}
\alias{assignTypes}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Assign genotypes to the query sequence
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
assignTypes("pathToRef", "pathToQuery", model = "p-distance", threshold = 0.105)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
 \item{pathToRef}{
  %%     ~~Describe \code{pathToRef} here~~
  pathToRef is path to the reference sequence
  }
  \item{pathToQuery}{
  %%     ~~Describe \code{pathToQuery} here~~
  pathToRef is path to the target sequence
  }
  \item{model}{
  %%     ~~Describe \code{model} here~~
  It is used to select the model used to assign genotypes to the query sequence. Default is set to p-distance.
  }
  \item{threshold}{
  %%     ~~Describe \code{threshhold} here~~
  The default threshhold is set to 0.105
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
  
  # run command
  fastaD <- readFasta(test)
  assignTypes(fastaD, "p-distance")
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
