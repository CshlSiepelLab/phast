\name{summary.msa}
\alias{summary.msa}
\title{
  MSA Summary
}
\description{
  Prints a short desription of an MSA (multiple sequence alignment) object.
}
\usage{
summary.msa(msa, ...)
}
\arguments{
  \item{msa}{
    MSA object.
}
  \item{\dots}{
    Additional arguments sent to \code{print.msa}.
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
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
m <- msa.new(names=c("human", "mouse", "rat"),
             seqs=c("ACGTAT", "AGGTAA", "AGGTAG"))
summary(m)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{print.msa}
\keyword{summary}% __ONLY ONE__ keyword per line
