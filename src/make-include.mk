###########################################################################
# this file defines variables used by all Makefiles
###########################################################################

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
CFLAGS = -g -fno-inline -Wall -pthread
# for best performance
#CFLAGS = -O3  -pthread
# some other options
#CFLAGS = -mcpu=opteron -O3
#CFLAGS = -mcpu=pentiumpro -O3 

CFLAGS += -I${INC} -DPHAST_VERSION=\"$(shell cat ${PHAST}/version)\" -DPHAST_HOME=\"${PHAST}\"
LIBPATH = -L${LIB} 

# uncomment these lines for profiling (add -g for line-by-line
# profiling and -a for monitoring of basic blocks)
#CFLAGS += -pg

# this flag tells certain routines to dump internal, debugging output.
# Don't uncomment unless you know what you're doing.
#CFLAGS += -DDEBUG

# uncomment this line to build the RPHAST shared library, which allows
# access to PHAST from the R programming environment.  Requires an
# up-to-date installation of R.
#RPHAST = T

ifdef RPHAST
CFLAGS += -fPIC
endif

# The next section is concerned with the LAPACK linear algebra
# package, which is used by PHAST for matrix diagonalization and
# matrix inversion.  You have two options: (1) If you are running Mac
# OS version 10.3 (Panther) or later, you can use the LAPACK libraries
# that are pre-installed as part of the vecLib framework; or (2) you
# can separately install the CLAPACK package and use its libraries
# (see README.txt for details).  You can also bypass LAPACK
# altogether, but in this case several key programs (including
# phastCons, exoniphy, and phyloFit) will not be usable.

# vecLib on Mac OS X; uncomment to use
#VECLIB = T

# separately installed CLAPACK; uncomment CLAPACKPATH definition and
# set appropriately to use
CLAPACKPATH = /usr/local/software/CLAPACK
# platform-specific suffix used for CLAPACK libraries; use the same
# value as in CLAPACK's "make.inc" file 
PLAT = _x86
# F2C libraries used by CLAPACK; most users won't need to edit
F2CPATH = ${CLAPACKPATH}/F2CLIBS


# if neither VECLIB nor CLAPACKPATH is defined, then LAPACK will be
# bypassed altogether


# Most users shouldn't edit the lines below (but see note about older
# versions of CLAPACK)

# vecLib
ifdef VECLIB
CFLAGS += -DVECLIB
LIBS = -lphast -framework vecLib -lc -lm

# CLAPACK
else
ifdef CLAPACKPATH
CFLAGS += -I${CLAPACKPATH}/INCLUDE -I${F2CPATH}
LIBS = -lphast -llapack -ltmg -lblaswr -lc -lf2c -lm
# IMPORTANT: use the following two lines instead for versions of CLAPACK
# older than 3.1.1
#CFLAGS += -I${CLAPACKPATH} -I${F2CPATH}
#LIBS = -lphast -llapack -ltmg -lblaswr -lc -lF77 -lI77 -lm
LIBPATH += -L${F2CPATH} 

# bypass
else
CFLAGS += -DSKIP_LAPACK
LIBS = -lphast -lc -lm
endif
endif

