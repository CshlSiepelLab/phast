---------------------------------------------------------------------------
            PHAST: PHYLOGENETIC ANALYSIS WITH SPACE/TIME MODELS
---------------------------------------------------------------------------

INTRODUCTION

This README describes a new software package for the analysis of multiple,
aligned biological sequences using combined phylogenetic and hidden Markov
models.  

The theory underlying PHAST is described in:

    A. Siepel and D. Haussler.  Combining phylogenetic and hidden Markov
      models in biosequence analysis.  RECOMB, 2003.

    A. Siepel and D. Haussler.  Phylogenetic estimation of
      context-dependent substitution rates by maximum likelihood.
      Mol. Biol. Evol., 21:468-488, 2004.

    A. Siepel and D. Haussler.  Computational identification of
      evolutionarily conserved exons.  RECOMB, 2004.

    A. Siepel and D. Haussler.  Phylogenetic hidden Markov models.  In
      R. Nielsen, Statistical Methods in Molecular Evolution.  In press.


DEVELOPMENT AND AVAILABILITY

PHAST was developed at the University of California, Santa Cruz, by
Adam Siepel, working under the direction of David Haussler (Computer
Science Department and Center for Biomolecular Science and Engineering).
All code has been developed in C under Linux, using GCC.  It has so far
been compiled and run only under Linux, but should be very portable.  Note
that a couple of auxiliary software packages are required (see below).

The software is freely available for research purposes.


CONTENTS

PHAST consists of a fairly large library of reusable subroutines and
a small collection of executable programs that use the library.
The library is composed of six packages, as follows:

    PACKAGE             CONTENTS
    ------------------------------------------------------------------------
    base                Basic data structures and supporting routines not
                        specific to phylogeny or HMMs (e.g., lists, strings,
                        matrices).

    msa                 Support for multiple sequence alignments

    feature             Code to map between sequence annotations and
                        labelings of alignment sites. 

    hmm                 HMM routines, including implementations of the
                        Viterbi and forward/backward algorithms.

    phylo               Support for phylogenetic models.  A variety of DNA
                        substitution models are supported.  Parameter
                        estimation can be accomplished by EM as well as by
                        standard quasi-Newton methods.

    phylo_hmm           Support for combined phylogenetic and hidden Markov
                        models.


(say something about executables)

INSTALLATION

Download and unpack the compressed archive file "phast.v*.tgz", using commands
such as "tar xfz phast.tgz" (with GNU tar) or "gunzip phast.tgz" and "tar
xf phast.tar".  A directory called "phast" will be created, containing this
README, as well as directories for source code (src), header files
(include), and documentation (doc).

PHAST currently depends on two auxiliary software packages that are not
available by default on most systems: the GNU Scientific Library (GSL) and
the CLAPACK linear algebra package.  You must download and install these
packages if you do not already have them (they are available at
http://www.gnu.org/software/gsl/ and http://www.netlib.org/clapack/,
respectively).  Edit the file "make-include.mk" in the "src" directory, so
that the variables "GSLPATH" and "CLAPACKPATH" are defined appropriately
(see the examples in the file).  You'll also need to either set a PHAST
environment variable to point to the root directory of your installation or
edit the variable of the same name in make-include.mk.  Note that PHAST
also depends on the "F2C" (Fortran to C) package, which is used by CLAPACK,
and on an implementation of the "BLAS" (Basic Linear Algebra Subroutines).
By default it uses the versions of these that come with CLAPACK.
Currently, linear algebra computations appear not to be a performance
bottleneck, so it's probably not worth going to much trouble to switch to a
highly optimized BLAS implementation.

The software also requires GNU Make, some standard UNIX tools such as "ar" and
"ln", and the GNU Regex package.  All of these should be available on most
UNIX systems, or via the Cygwin toolkit (or similar) for Windows.

Once the GSL and CLAPACK packages have been installed and "make-include.mk"
has been editted, you can simply "cd" to the "src" directory and type
"make".  The library and executables should compile cleanly.  If not,
please report problems to acs@soe.ucsc.edu.

The libraries and executables will be placed in directories called "lib"
and "bin", respectively, at the top level of the "phast" tree. 

