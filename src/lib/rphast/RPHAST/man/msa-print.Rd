\name{print.msa}
\alias{print.msa}
\title{
  Printing MSA objects.
  %%  ~~function to do ... ~~
}
\description{
  Prints an msa (multiple sequence alignment) object.
}
\usage{
print.msa(msa, ..., printSeq = FALSE, format = "FASTA", prettyPrint = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{msa}{
    An object of class 'msa'.
  }
  \item{\dots}{optional arguments to \code{print} methods.}
  \item{printSeq}{
    If TRUE, does not supress the alignment.
  }
  \item{format}{
    If printSeq is TRUE, specifies the format to print the alignment.
    See \code{\link{msa.validFormat}}.
  }
  \item{prettyPrint}{
    If printSeq is TRUE, prettyPrint=TRUE will cause '.' to be printed
    whenever character matches corresponding character in first
    sequence.  Ignored if format is SS.
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
m
print(m, printSeq=TRUE)
print(m, printSeq=TRUE, format="SS")
print(m, printSeq=TRUE, format="PHYLIP", prettyPrint=TRUE)

  }
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{print}
\keyword{msa}% __ONLY ONE__ keyword per line
