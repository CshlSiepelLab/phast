---------------------------------------------------------------------------
            PHAST: PHYLOGENETIC ANALYSIS WITH SPACE/TIME MODELS
---------------------------------------------------------------------------

This README describes a software package for the analysis of multiple,
aligned biological sequences using combined phylogenetic and hidden
Markov models.

The theory underlying PHAST is described in:

    A. Siepel and D. Haussler.  Combining phylogenetic and hidden Markov
      models in biosequence analysis.  RECOMB, 2003.

    A. Siepel and D. Haussler.  Phylogenetic estimation of
      context-dependent substitution rates by maximum likelihood.
      Mol. Biol. Evol., 21:468-488, 2004.

    A. Siepel and D. Haussler.  Computational identification of
      evolutionarily conserved exons.  RECOMB, 2004.

    A. Siepel and D. Haussler.  Phylogenetic hidden Markov models.  In
      R. Nielsen, Statistical Methods in Molecular Evolution,
      Springer, NY, 2005.


		     DEVELOPMENT AND AVAILABILITY

PHAST was developed at the University of California, Santa Cruz, by
Adam Siepel, working under the direction of David Haussler (Computer
Science Department and Center for Biomolecular Science and
Engineering).  The package is written in ANSI C, and has been compiled
and run on the Linux, SunOS, and Mac OS X operating systems; it should
be easy to port it to other UNIX-based systems.  Note that it
currently depends on the CLAPACK linear algebra package (see below).

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



		      INSTALLATION INSTRUCTIONS

1. Download and unpack the compressed archive file "phast.v*.tgz",
   using commands such as "tar xfz phast.tgz" (with GNU tar) or
   "gunzip phast.tgz" and "tar xf phast.tar".  A directory called
   "phast" will be created, containing this README and directories for
   source code (src), header files (include), and documentation (doc).

2. Install the CLAPACK linear algebra package, which is used by PHAST
   for certain matrix manipulations (viz. diagonalization).

    a. Download clapack.tgz from http://www.netlib.org/clapack, unpack
       in your directory of choice (e.g., "tar xfz clapack.tgz"), and cd
       to the CLAPACK directory

    b. Edit the following variables in make.inc.

          PLAT: Set equal to an appropriate string (any will do) for
            your platform, e.g., "PLAT = _x86" or "PLAT = _MACOSX"

          CFLAGS: Set appropriate compiler options, e.g.,
	    "CFLAGS = -mcpu=pentiumpro -O3" for x86 or "CFLAGS =
	    -mcpu=powerpc -O3 -Wno-long-double" for Mac OS X.

          LOADOPTS: "LOADOPTS = $(CFLAGS)" will do.
    
    c. Follow the quick-start instructions in README.install.  You
       should just be able to type "make f2clib", then "make blaslib",
       then "make".

3. Edit the file "make-include.mk" in phast/src.  Specifically:

    a. Set the PHAST variable equal to the root directory of your
       installation (the "phast" directory)

    b. Set CFLAGS appropriately.  You may just be able to uncomment
       one of the lines in the file.

    c. Set the CLAPACKPATH variable equal to the root directory of
       your CLAPACK installation.

    d. Set the ARCH variable to the string you used for PLAT with
       CLAPACK (e.g., "_x86" or "_MACOSX").  

    e. With MAC OS X, also comment out the "LFLAGS = -static" option.

4. Type "cd src" and "make".  The package should compile cleanly.  If
   you encounter problems compiling, please report them to
   acs@soe.ucsc.edu.  I'll do my best to help you work around them and
   to avoid similar problems in the future.


				NOTES

    - PHAST also depends on the "F2C" (Fortran to C) package and on an
      implementation of "BLAS" (Basic Linear Algebra Subroutines),
      both of which are used by CLAPACK.  By default it uses the
      versions of these that come with CLAPACK.  The default BLAS
      implementation seems to be fine for normal usage.

    - The software requires GNU Make, some standard UNIX tools (e.g.,
      sed, ar, and ln), and a getopt implementation that supports long
      options (e.g., GNU getopt or BSD getopt).  These should be
      available on most UNIX systems, on Mac OS X, and via the Cygwin
      toolkit for Windows.

    - It's possible to compile the software without CLAPACK by
      uncommenting the SKIP_CLAPACK line in make-include.mk.  In this
      case, some programs will be usable, but programs that require
      matrix diagonalization will abort at the critical point of
      calling a CLAPACK routine.

