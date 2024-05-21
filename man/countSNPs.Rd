\name{countSNPs}
\alias{countSNPs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Counts single nucleotide polymorphisms.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Counts single nucleotide polymorphisms across input sequences and produce an output matrix.
}
\usage{
countSNPs(inputSequencesPath)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{inputSequencesPath}{
%%     ~~Describe \code{inputSequencesPath} here~~
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
  
  # run
  fastaData <- readFasta(test)
  countSNPs(fastaData, gapDeletion = TRUE))
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
