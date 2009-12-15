phast.init <- function(){
   dyn.load("/home/melissa/phast/lib/librphast.so")
   return()
}

msa.free <- function(extMsaPtr) {
  print("msa.free")
  .Call("rph_msa_free", extMsaPtr)
}

msa.new <- function(nseqs=0, seqlen=0, names=NULL, seqs=NULL, alphabet=NULL) {
  msa <- list()
  class(msa) <- "msa"
  if (! is.null(seqs)) {
    if (seqlen == 0) {
      tempSeqLens <- unique(sapply(seqs, nchar))
      if (length(tempSeqLens) > 1) 
        stop("sequences not all equal length")
      seqlen<-tempSeqLens
    }
    if (nseqs == 0) 
      nseqs<-length(seqs)
    if (!is.null(names) & length(names) != nseqs)
      stop("number of names does not match number of sequences")
  }
  
  msa$nseqs  <- nseqs
  msa$seqlen <- seqlen
  msa$names <- names
  msa$alphabet <- alphabet

  msa$externalPtr <- .Call("rph_msa_new",
                           seqsP=seqs,
                           namesP=names,
                           nseqsP=nseqs,
                           lengthP=seqlen,
                           alphabetP=alphabet)
  reg.finalizer(msa$externalPtr, msa.free)
  msa
}

msa.validFormat <- function(format) {
  isValid <- .Call("rph_msa_valid_fmt", format);
  isValid
}


print.msa <- function(msa, printSeq=FALSE, format="FASTA", prettyPrint=FALSE) {
  cat(paste("msa object with", msa$nseqs, "sequences and", msa$seqlen, "columns"))
  cat("\n")
  if (!is.null(msa$names)) {
    cat("$names\n")
    print(msa$names)
  }
  if (!is.null(msa$alphabet)) {
    cat("$alphabet")
    print(msa$alphabet)
  }
  if (printSeq) {
    if (! msa.validFormat(format)) {
      warning(paste("invalid MSA FORMAT \"", format, "\"", sep=""))
    } else {
      x<-.Call("rph_msa_printSeq",
               msaP=msa$externalPtr,
               formatP=format,
               prettyPrintP=(prettyPrint==TRUE))
    }
  }
}

# same as print.msa but no option to print sequence
summary.msa <- function(msa) {
  print.msa(msa)
}


########################################
#testing code
msa.testing <- function() {
seqs <- c("ACGT", "ACGT", "ACCC")
names <- c("human", "chimp", "rat")
nseqs <- 3
length <- 4
alphabet <- "ACGT"

dyn.load("/home/melissa/phast/lib/librphast.so")
m <- msa.new(3, 4, names, seqs, NULL)
rm(m)
gc()

m <- msa.new(3, 4, names, seqs, NULL)
print(m)
print(m, printSeq=TRUE)
print(m, printSeq=TRUE, format="PHYLIP")
}
