\name{plotAA}
\alias{plotAA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Plots amino acid substitutions with a specified sequence as the reference. The input is an amino acid fasta file (translated DNA sequences). To specify the reference sequence, move it to the bottom of the alignment. Changes are coloured by the class of amino acid: Red = Positively charged, Blue = Negatively charged, Green = Polar, Yellow = Non-polar.
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
