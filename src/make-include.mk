###########################################################################
# this file defines variables used by all Makefiles
###########################################################################

# currently set up for linux and GCC and Pentium III/IV processors,
# with options for Mac OS X.  Adjustments may be needed for other
# systems.

# set below to point to top-level directory of PHAST installation
ifndef PHAST
PHAST=${HOME}/phast
endif
# (if you prefer, you can set the environment variable PHAST instead)

# specify alternative compiler or utilities if necessary
CC = gcc
AR = ar
LN = ln

LIB = ${PHAST}/lib
INC = ${PHAST}/include
BIN = ${PHAST}/bin

TARGETLIB = ${LIB}/libphast.a

# set compiler options; uncomment one of the lines below or define
# an appropriate alternative

# for debugging
#CFLAGS = -g -fno-inline -Wall 
# for best performance (pentiumpro)
CFLAGS = -mcpu=pentiumpro -O3 
# use this instead for Mac OS X
#CFLAGS = -mcpu=powerpc -O3
# possible x86-64 options (kolossus at UC Santa Cruz)
#CFLAGS = -mcpu=opteron -O3

CFLAGS += -I${INC} -DPHAST_VERSION=\"$(shell cat ${PHAST}/version)\" -DPHAST_HOME=\"${PHAST}\"
LIBPATH = -L${LIB} 

# static linking ends up being simplest in our environment; comment
# this line out to link dynamically (comment out for Mac)
LFLAGS += -static

# comment these lines out for profiling (add -g for line-by-line
# profiling and -a for monitoring of basic blocks)
CFLAGS += -pg
LFLAGS += -pg

# Uncomment this line to compile without CLAPACK.  If SKIP_CLAPACK is
# defined, PHAST will compile, but programs that require matrix
# diagonalization will not be usable.   
#SKIP_CLAPACK = T

# root directory for CLAPACK; necessary unless SKIP_CLAPACK is defined
CLAPACKPATH = /projects/compbio/usr/acs/CLAPACK

# location of F2C libraries used by CLAPACK; most people won't need to
# edit this
F2CPATH = ${CLAPACKPATH}/F2CLIBS

# platform-specific suffix used for CLAPACK libraries; use the same
# value used in the "make.inc" file when compiling CLAPACK
PLAT = _x86

ifndef SKIP_CLAPACK
CFLAGS += -I${CLAPACKPATH} -I${F2CPATH}
LIBPATH += -L${F2CPATH} 
LIBS = -lphast -llapack -ltmg -lblaswr -lc -lF77 -lI77 -lm
else
CFLAGS += -DSKIP_CLAPACK
LIBS = -lphast -lm
endif

# this flag tells certain routines to dump internal, debugging output.
# Don't uncomment unless you know what you're doing.
#CFLAGS += -DDEBUG

