\name{msa.validFormat}
\alias{msa.validFormat}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Check an MSA format string
%%  ~~function to do ... ~~
}
\description{
  Returns TRUE if the argument is a valid format string for a multiple
  sequence alignment (MSA) object.
}
\usage{
msa.validFormat(format)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{format}{
    the character string to test.  Can also be a vector of character
    strings.
%%     ~~Describe \code{format} here~~
}
}
\details{
  Valid formats include "FASTA", "PHYLIP", "SS" (Sufficient statistics
  format used by PHAST), "MPM" (format used by MultiPipMaker), "LAV"
  (used by blastz), or "MAF" (Multiple Alignment Format used by MULTIZ
  and TBA).
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  Returns a logical vector.  Each element is TRUE if the corresponding format matches a valid format, otherwise FALSE.
  
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
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
msa.validFormat(c("MAF", "SS", "PHYLIP", "MPM", "LAV", "FASTA", "BAD_FORMAT_STRING"))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
