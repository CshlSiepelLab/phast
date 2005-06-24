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
#CFLAGS = -g -fno-inline -Wall -DGCC -Wno-long-double
# for best performance (pentiumpro)
CFLAGS = -mcpu=pentiumpro -O3 -DGCC
# use this instead for Mac OS X
#CFLAGS = -mcpu=powerpc -O3 -DGCC -Wno-long-double
# possible x86-64 options (kolossus at UC Santa Cruz)
#CFLAGS = -mcpu=opteron -O3 -DGCC
# possible profiling options
#CFLAGS = -mcpu=pentiumpro -O3 -DGCC -pg -fprofile-arcs
# NOTE: add -g if profiling line-by-line, and -a if monitoring basic blocks

CFLAGS += -DPHAST_VERSION=\"$(shell cat ${PHAST}/version)\" -DPHAST_HOME=\"${PHAST}\"

# static linking ends up being simplest in our environment; comment
# this line out to link dynamically (comment out for Mac)
LFLAGS = -static

# Uncomment this line to compile without CLAPACK.  If SKIP_CLAPACK is
# defined, PHAST will compile, but programs that require matrix
# diagonalization will not be usable.   
#SKIP_CLAPACK = T

# root directory for CLAPACK; necessary unless SKIP_CLAPACK is defined
CLAPACKPATH = /projects/compbio/usr/acs/CLAPACK

# location of F2C libraries used by CLAPACK; you shouldn't need to
# edit this unless you're using a different version of F2C
F2CPATH = ${CLAPACKPATH}/F2CLIBS

# platform-specific suffix used for CLAPACK libraries; use the same
# value used in the "make.inc" file when compiling CLAPACK
PLAT = _x86

CFLAGS += -I${INC}
LIBPATH = -L${LIB} 

ifndef SKIP_CLAPACK
CFLAGS += -I${CLAPACKPATH}
ifdef F2CPATH
CFLAGS += -I${F2CPATH}
LIBPATH += -L${F2CPATH} 
endif
LIBS = -lphast -llapack -ltmg -lblaswr -lc -lF77 -lI77 -lm
else
CFLAGS += -DSKIP_CLAPACK
LIBS = -lphast -lm
endif

# this flag tells certain routines to dump internal, debugging output.
# Don't uncomment unless you know what you're doing.
#CFLAGS += -DDEBUG

