###########################################################################
# this file defines variables used by all Makefiles
###########################################################################

# currently set up for linux and GCC and Pentium III/IV processors.
# Adjustments will be needed for other systems.

# set either environment variable PHAST or path below to top-level directory
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
#CFLAGS = -g -fno-inline -Wall -DGCC
# for best performance
CFLAGS = -mcpu=pentiumpro -O3 -DGCC -DGSL_RANGE_CHECK_OFF
# for profiling
#CFLAGS = -mcpu=pentiumpro -O3 -DGCC -pg -fprofile-arcs -DGSL_RANGE_CHECK_OFF
# NOTE: add -g if profiling line-by-line, and -a if monitoring basic blocks

CFLAGS += -DPHAST_VERSION=\"$(shell cat ${PHAST}/version)\" -DPHAST_HOME=\"${PHAST}\"

# static linking ends up being simplest in our environment; comment
# this line out to link dynamically
LFLAGS = -static

# define if GSL files not in standard locations (such as /usr/include/gsl
# or /usr/lib)
GSLPATH = /projects/compbio/usr/acs/gsl-1.4-P3

# must be defined (see Makefile in lib subdir)
CLAPACKPATH = /projects/compbio/usr/acs/CLAPACK

# define if F2C files not in standard locations
F2CPATH = ${CLAPACKPATH}/F2CLIBS

# architecture-specific suffix used for CLAPACK libraries; acceptable
# values include LINUX, SUN4SOL2, HPPA, ALPHA, OCTANE, RS6K
ARCH = LINUX

CFLAGS += -I${INC} -I${CLAPACKPATH}
LIBPATH = -L${LIB} 
ifdef GSLPATH
	CFLAGS += -I${GSLPATH}/include
	LIBPATH += -L${GSLPATH}/lib
endif
ifdef F2CPATH
	CFLAGS += -I${F2CPATH}
	LIBPATH += -L${F2CPATH} 
endif

LIBS =  -lphast -lgsl -llapack -ltmg -lblaswr -lc \
	-lF77 -lI77 -lgslcblas -lm

# NOTE: assuming version of BLAS that comes with CLAPACK; might want
# to allow alternatives.  Also, note that gslcblas is required because
# of use of BLAS via GSL wrappers

# this flag tells certain routines to dump internal, debugging output.  
#CFLAGS += -DDEBUG

