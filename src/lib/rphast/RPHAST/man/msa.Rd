\name{msa.new}
\alias{msa.new}
\alias{msa}
\title{MSA ojbects}
\description{
  Creates a new MSA (multiple sequence alignment) object.
}
\usage{
msa.new(nseqs = 0, seqlen = 0, names = NULL, seqs = NULL, alphabet = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{nseqs}{
    the number of sequences (ie, samples or rows) in the alignment.
}
\item{seqlen}{
  the length of the alignment (usually the number of base pairs).
}
  \item{names}{
    a character vector of length nseqs with the name of each sample,
    in the same order as given in the alignment.

}
  \item{seqs}{
    a character vector of length nseqs, each string of length seqlen,
    containing the alignment.
}
  \item{alphabet}{
    a character string containing all valid characters that may appear
  in the alignment.
}
}
\details{
  If seqs is not NULL, nseqs and seqlen will be inferred from seqs, and
  a warning will be given if they do not match the arguments provided.
  The alphabet will not be inferred from seqs unless it is required in a
  downstream operation.
}
\value{
  Returns a new msa object.
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

\section{Warning }{ Unlike most R objects, msa objects do not "copy by
  value".  Due to the potentially large size of alignments, the
  alignment is not copied unless msa.copy() is called explicitly.  An
  assignment msa2 <- msa1 will result in both MSA objects pointing to
  the same alignment in memory.}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  \code{\link{msa.copy}}, \code{\link{print.msa}}, \code{\link{write.msa}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{

m <- msa.new(nseqs=3, seqlen=6, names=c("human", "mouse", "rat"))
m
print(m, printSeq=TRUE)

m <- msa.new(names=c("human", "mouse", "rat"),
             seqs=c("ACGTAT", "AGGTAA", "AGGTAG"))
m
print(m, printSeq=TRUE)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
