---------------------------------------------------------------------------
PHAST: PHYLOGENETIC ANALYSIS WITH SPACE/TIME MODELS
---------------------------------------------------------------------------

QUICK START - INSTALLING PHAST

You can either download the Phast binaries or compile Phast from source.
Binaries are generally the easiest way to get up and running with Phast 
and are suggested for new users.  

-Installing using Binaries
  Phast binaries for Windows, MacOSX, and Linux can be downloaded from the 
    main PHAST page: http://compgen.bscb.cornell.edu/phast

-Compiling from Source	
  MacOSX
    1. Download a copy of Phast from http://compgen.bscb.cornell.edu/phast/ 
	and extract the file phast*.tgz using 'tar -xvzf phast*.tgz'
    2. Change directory to 'phast/src/' and run the command 'make'
    3. The Phast binaries should be located in the 'phast/bin/' directory.

  Linux
    Part 1 - Installing Clapack - (If you already have Clapack installed, skip to Part 2)
    1. Download Clapack from the following URL http://www.netlib.org/clapack/clapack.tgz
    2. Unzip clapack.tgz with the command 'tar -xvzf clapack.tgz'
    3. Go into the newly created Clapack directory (i.e. 'cd CLAPACK-3.2.1') 
	and type 'cp make.inc.example make.inc && make f2clib && make blaslib && make lib'
       Note: Building Clapack can take several minutes depending on your system
    
    Part 2 - Installing Phast
    4. Download a copy of Phast from http://compgen.bscb.cornell.edu/phast/ 
	and extract the contents of phast*.tgz using 'tar -xvzf phast*.tgz'
    5. Change directory to 'phast/src/' and run 'make CLAPACKPATH=/usr/local/software/clapack'
	replacing '/usr/local/software/clapack' with the path of your
	Clapack install (e.g., CLAPACKPATH=/home/username/CLAPACK-3.2.1)
    6. The Phast binaries should be created in the '../bin/' directory
    
  Windows 
    PHAST can be compiled under Windows using the Cygwin linux-like
    environment, but some users have reported difficulties in making this
    work.  We recommend using the provided binaries for Windows, unless you
    have a good reason to compile the package from source.

The Phast package should compile cleanly in most standard linux or
linux-like environments (including MacOS).  If you encounter problems
compiling, please report them to phast-help-l@cornell.edu.  We'll do our
best to help you work around them and to avoid similar problems in the
future.


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

    - The most recent source code of Phast can be obtained from our public
      subversion server.  If you are set up to use subversion, you may want
      to check out the latest version before submitting a bug report.
      To download the latest code simply run the command as follows
      'svn co http://compgen.bscb.cornell.edu/svnrepo/phast/trunk/phast'.
      The source code will be saved in a folder called 'phast' in the 
      directory where the command was run. You can substitute SVN code for
      the *.tgz formated code download from our site and build by following
      the above instructions.


ACKNOWLEDGEMENTS

PHAST makes use of the CLAPACK linear algebra library
(http://www.netlib.org/clapack/) and the PCRE regular expression library
(http://www.pcre.org).  We thank the authors of these packages for making
them freely available to the community.
