###########################################################################
# this file defines variables used by all Makefiles
###########################################################################

# currently set up for linux and GCC and Pentium III/IV processors,
# with options for Mac OS X.  Adjustments will be needed for other
# systems.

# set either below or environment variable PHAST to point to top-level directory
ifndef PHAST
PHAST=${HOME}/phast
endif

CC = gcc
AR = ar
LN = ln
LIB = ${PHAST}/lib
INC = ${PHAST}/include
BIN = ${PHAST}/bin

TARGETLIB = ${LIB}/libphast.a

# for debugging
CFLAGS = -g -fno-inline -Wall -DGCC -Wno-long-double
# for best performance
#CFLAGS = -mcpu=pentiumpro -O3 -DGCC -DGSL_RANGE_CHECK_OFF
# use this instead for Mac OS X
#CFLAGS = -mcpu=powerpc -O3 -DGCC -DGSL_RANGE_CHECK_OFF -Wno-long-double
# possible x86-64 options (kolossus at UC Santa Cruz)
#CFLAGS = -mcpu=opteron -O3 -DGCC -DGSL_RANGE_CHECK_OFF
# possible profiling options
#CFLAGS = -mcpu=pentiumpro -O3 -DGCC -pg -fprofile-arcs -DGSL_RANGE_CHECK_OFF
# NOTE: add -g if profiling line-by-line, and -a if monitoring basic blocks

CFLAGS += -DPHAST_VERSION=\"$(shell cat ${PHAST}/version)\" -DPHAST_HOME=\"${PHAST}\"

# static linking ends up being simplest in our environment; comment
# this line out to link dynamically (comment out for Mac)
#LFLAGS = -static

# define if GSL files not in standard locations (such as /usr/include/gsl
# or /usr/lib)
#GSLPATH = /projects/compbio/usr/acs/gsl-1.4-P3
#GSLPATH = /Users/acs/UTILITIES/gsl-1.6
# above to be removed

# Uncomment this line to compile without CLAPACK.  If SKIP_CLAPACK is
# defined, PHAST will compile, but programs that require matrix
# diagonalization will not be usable.   
#SKIP_CLAPACK = T

# must be defined unless SKIP_CLAPACK (see Makefile in lib subdir)
#CLAPACKPATH = /projects/compbio/usr/acs/CLAPACK
CLAPACKPATH = /Users/acs/UTILITIES/CLAPACK

# define if no SKIP_CLAPACK and F2C files not in standard locations
F2CPATH = ${CLAPACKPATH}/F2CLIBS

# architecture-specific suffix used for CLAPACK libraries; acceptable
# values include LINUX, SUN4SOL2, HPPA, ALPHA, OCTANE, RS6K, MACOSX
#ARCH=LINUX
ARCH=MACOSX

CFLAGS += -I${INC}
LIBPATH = -L${LIB} 

# to be removed
#ifdef GSLPATH
#CFLAGS += -I${GSLPATH}/include
#LIBPATH += -L${GSLPATH}/lib
#endif

ifndef SKIP_CLAPACK
CFLAGS += -I${CLAPACKPATH}
ifdef F2CPATH
CFLAGS += -I${F2CPATH}
LIBPATH += -L${F2CPATH} 
endif
LIBS =  -lphast -llapack -ltmg -lblaswr -lc \
	-lF77 -lI77 -lm
#LIBS =  -lphast -lgsl -llapack -ltmg -lblaswr -lc \
#	-lF77 -lI77 -lgslcblas -lm
else
CFLAGS += -DSKIP_CLAPACK
LIBS =  -lphast -lc -lm
endif

# remove gsl and gslcblas above; is -lc necessary?

# NOTE: currently assuming version of BLAS that comes with CLAPACK;
# might want to allow alternatives.

# this flag tells certain routines to dump internal, debugging output.  
#CFLAGS += -DDEBUG

