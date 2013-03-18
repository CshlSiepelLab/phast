###########################################################################
# this file defines variables used by all Makefiles
###########################################################################

# If the user did not specify a Operating System to target, determine what OS this system is using 
ifndef TARGETOS
  TARGETOS := $(shell uname -s)
endif
# (if you prefer, you can target a specific OS instead by setting the environment variable TARGETOS instead)

# Points to top-level directory of PHAST installation
ifndef PHAST
  PHAST=${PWD}/..
endif

#If the user does not specify PHAST_HOME then set it to the old deault of ${PHAST}
ifndef PHAST_HOME
  PHAST_HOME=${PHAST}
endif
# (if you prefer, you can set the environment variable PHAST instead)

# specify alternative compiler or utilities if necessary
ifeq ($(TARGETOS), Windows)
  CC = /usr/bin/i586-mingw32msvc-gcc
  AR = /usr/bin/i586-mingw32msvc-ar
else
  ifeq ($(TARGETOS), LSB)
    CC = lsbcc -fno-stack-protector    
  else
    CC = gcc
  endif
  AR = ar
endif
LN = ln

LIB = ${PHAST}/lib
INC = ${PHAST}/include
BIN = ${PHAST}/bin

TARGETLIB = ${LIB}/libphast.a

# set compiler options; uncomment one of the lines below or define
# an appropriate alternative
ifneq ($(TARGETOS), Windows)
 #for debugging
# CFLAGS = -g -fno-inline -Wall -DPHAST_DEBUG
 # for best performance
 CFLAGS = -O3 
 # some other options
 #CFLAGS = -mcpu=opteron -O3
 #CFLAGS = -mcpu=pentiumpro -O3 
else
  CFLAGS = -O3
endif

PHAST_VERSION=\"$(shell cat ${PHAST}/version)\"
CFLAGS += -I${INC} -DPHAST_VERSION=${PHAST_VERSION} -DPHAST_HOME=\"${PHAST_HOME}\" -I${PHAST}/src/lib/pcre -fno-strict-aliasing
LIBPATH = -L${LIB} 

# uncomment these lines for profiling (add -g for line-by-line
# profiling and -a for monitoring of basic blocks)
#CFLAGS += -pg

# this flag tells certain routines to dump internal, debugging output.
# Don't uncomment unless you know what you're doing.
#CFLAGS += -DDEBUG

# ignore the section below if installing RPHAST
ifndef RPHAST

# the following line should be uncommented during phast development
# to make sure RPHAST C files stay up-to-date.
# Not necessary for installing phast or RPHAST.
#RPHAST = T
ifdef RPHAST
RDIR=/usr/share/R/include
CFLAGS += -fPIC -I${RDIR}
endif
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
ifeq ($(TARGETOS), Darwin)
  VECLIB = T
endif

# separately installed CLAPACK; uncomment CLAPACKPATH definition and
# set appropriately to use, or define CLAPACKPATH when you run 'make'
ifndef VECLIB
  # platform-specific suffix used for CLAPACK libraries; use the same
  #value as in CLAPACK's "make.inc" file
ifneq ($(TARGETOS), Windows)
  ifndef CLAPACKPATH
    CHECKFILE = $(shell if [ -d /usr/local/software/clapack ]; then echo "true"; fi) 
    ifeq ($(CHECKFILE),true )
      CLAPACKPATH = /usr/local/software/clapack
    endif 
  endif 
    ifndef CLAPACKPATH
      CLAPACKPATH = /usr/local/software/clapack
    endif 
    #Automatically detects PLAT type by looking in CLAPACKPATH for blas*.a and extracts the * part
    PLAT = $(shell find ${CLAPACKPATH}/ -name '*.a' -exec expr match {} '.*blas\(.*\).a' \; | tr -d "\n")
  else
    ifndef CLAPACKPATH
	    CLAPACKPATH = /usr/local/software/clapack-windows
    endif
    # PLAT is empty for windows builds
    PLAT =
  endif
  # F2C libraries used by CLAPACK; most users won't need to edit
  F2CPATH = ${CLAPACKPATH}/F2CLIBS
endif

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
ifneq ($(TARGETOS), Windows)
  CFLAGS += -I${CLAPACKPATH}/INCLUDE -I${F2CPATH}
  LIBS = -lphast -llapack -ltmg -lblaswr -lc -lf2c -lm
else
  CFLAGS += -I${CLAPACKPATH}/INCLUDE -I${F2CPATH} -DPCRE_STATIC
  LIBS = -lphast -lm  ${CLAPACKPATH}/liblapack.a ${CLAPACKPATH}/libf2c.a ${CLAPACKPATH}/libblas.a
endif
# IMPORTANT: use the following two lines instead for versions of CLAPACK
# older than 3.1.1
#CFLAGS += -I${CLAPACKPATH} -I${F2CPATH}
#LIBS = -lphast -llapack -ltmg -lblaswr -lc -lF77 -lI77 -lm
LIBPATH += -L${F2CPATH} 

# bypass
else
ifneq ($(TARGETOS), Windows)
  CFLAGS += -DSKIP_LAPACK
  LIBS = -lphast -lc -lm
else
  CFLAGS += -DSKIP_LAPACK -DPCRE_STATIC
  LIBS = -lphast -lm  
endif
endif
endif

