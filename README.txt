---------------------------------------------------------------------------
            PHAST: PHYLOGENETIC ANALYSIS WITH SPACE/TIME MODELS
---------------------------------------------------------------------------

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

2. If necessary, install the CLAPACK linear algebra package, which is used
   by PHAST for certain matrix manipulations (diagonalization and
   inversion).  Note that this step is NOT NECESSARY if using Mac OS
   version X.3 (Panther) or later (see comments in src/make-include.mk).  

   IMPORTANT NOTE: The organization of the CLAPACK package changed slightly
   in version 3.1.1 (released Feb 2008).  This version of PHAST assumes the
   new organization.  See instructions in src/make-include.mk if using an
   older version.

    a. Download clapack.tgz from http://www.netlib.org/clapack, unpack
       in your directory of choice (e.g., "tar xfz clapack.tgz"), and cd
       to the CLAPACK directory.

       Note: if the file clapack.h is not in your CLAPACK/INCLUDE/
       directory, you will need to download it separately from 
       http://www.netlib.org/clapack and put it there.  It seems that 
       some more recent versions of clapack do not include it in the
       tarball.

    b. Edit the following variables in make.inc.

          PLAT: Set equal to an appropriate suffix (any will do) for
            your platform, e.g., "PLAT = _x86" or "PLAT = _MACOSX"

          CFLAGS: Set appropriate compiler options, e.g.,
            "CFLAGS = -mcpu=pentiumpro -O3" for x86 or "CFLAGS =
            -mcpu=powerpc -O3" for Mac OS X.

          (IMPORTANT: for RPHAST, must include "-fPIC" among CFLAGS options
          when compiling CLAPACK)

          NOOPT: Leave undefined: "NOOPT = " 

          LOADOPTS: "LOADOPTS = $(CFLAGS)" will do.
    
    c. Follow the instructions in README.install, or, for a quick start,
       just type "make blaslib", then "make lib".  

       CLAPACK INSTALLATION NOTES: 

         (1) Under some conditions, a problem seems to occur on the Mac
             with the table of contents of the CLAPACK archive libraries.
             This can be remedied by running 'ranlib' manually on the
             affected files, as suggested by the error message.

         (2) A user has also reported CLAPACK installation errors under x86
             Linux similar to those documented at
             http://sourceware.redhat.com/ml/cygwin/2004-11/msg00666.html. He
             simply repeated the command 'cd CLAPACK/TESTING; make' several
             times (with different errors each run), until no errors
             remained.


3. Edit the file "make-include.mk" in phast/src.  Specifically:

    a. Set the PHAST variable equal to the root directory of your
       installation (the "phast" directory)

    b. Set CFLAGS appropriately.  You may just be able to uncomment
       one of the lines in the file.

    c. Following the directions in the file, set either the VECLIB variable
       (if Mac OS X.3 or later), or the CLAPACKPATH and PLAT variables, to
       allow linking to LAPACK.

4. Type "cd src" and "make".  The package should compile cleanly.  If
   you encounter problems compiling, please report them to
   phast-help-l@cornell.edu.  I'll do my best to help you work around
   them and to avoid similar problems in the future.


                                   NOTES

    - If CLAPACK is used, PHAST also depends on the "F2C" (Fortran to C)
      package and on an implementation of the "BLAS" (Basic Linear Algebra
      Subroutines).  By default it uses the versions of these that come
      with CLAPACK.  The default BLAS implementation seems to be fine for
      normal usage.

    - The software requires GNU Make, some standard UNIX tools (e.g.,
      sed, ar, and ln), and a getopt implementation that supports long
      options (e.g., GNU getopt or BSD getopt).  These should be
      available on most UNIX systems, on Mac OS X, and via the Cygwin
      toolkit for Windows.

    - It's possible to compile the software without LAPACK by commenting
      out both the VECLIB and CLAPACKPATH lines in src/make-include.mk.  In
      this case, some programs will be usable, but programs that require
      matrix diagonalization will abort at the critical point of calling a
      LAPACK routine.

