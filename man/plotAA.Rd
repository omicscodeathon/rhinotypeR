\name{plotAA}
\alias{plotAA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Visualise amino acid substitutions using a user-specified sequence as the reference
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
plotAA(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
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
  test <- system.file("extdata", "test.fasta", package = "rhinotypeR")
  
  # usage
  fastaData <- readFasta(fastaFile = test, desiredLength = 170)
  plotAA(fastaData, showLegend = T)
  
  # Legend
   # R = "red", H = "red", K = "red",      # Positively charged amino acid 
   # D = "blue", E = "blue",               # Negatively charged amino acid 
   # S = "green", T = "green", N = "green", Q = "green", # Polar amino acid 
   # A = "yellow", V = "yellow", I = "yellow", L = "yellow", M = "yellow", F = "yellow", 
   # W = "yellow", P = "yellow", G = "yellow", Y = "yellow", C = "yellow"  # Nonpolar amino acid 
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
