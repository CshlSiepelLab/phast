\name{msa.write}
\alias{msa.write}
\title{
  Writing MSA objects to files
}
\description{
  Writes a multiple sequence alignment (MSA) object to a file in one of
  several format.
}
\usage{
msa.write(msa, filename = NULL, format = "FASTA", prettyPrint = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{msa}{
    MSA object.
}
\item{filename}{
  File to write (will be overwritten).
}
\item{format}{
  Format.  Should be valid MSA format.  See \code{msa.validFormat}.  MAF
  format is not currently supported for output.
}
  \item{prettyPrint}{
    Use pretty-printing (replacing characters which match the first
    sequence to '.').  Ignored if format=="SS".
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
msa.write(m, "testSeq.ss", format="SS")
msa.write(m, "testSeq.fa", format="FASTA", prettyPrint=TRUE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{write}
\keyword{msa}% __ONLY ONE__ keyword per line
