#!/bin/bash -e
PHAST_VERSION=$1
CLAPACK=$2
PLAT=$3
RPHAST=$4

mkdir -p RPHAST/src/include
cp -ua `find ../ -name "*.c" | grep -v "RPHAST/src"` RPHAST/src
cp -ua ../../../include/ RPHAST/src

CSOURCES=`cd RPHAST/src; ls *.c`

# if development, make Makefile instead of Makevars
if [[ $RPHAST == "D" ]]; then
echo -e \
"include ../../../../make-include.mk
OBJECTS = \$(addsuffix .o, \$(basename \$(wildcard *.c)))
all: \$(OBJECTS)
\t\$(CC) -shared -o rphast.so \$(CFLAGS) \$(LIBPATH) \$(LFLAGS) -L/usr/lib64/R/lib -lR \$(OBJECTS) -I$CLAPACK/INCLUDE $CLAPACK/lapack_$PLAT.a $CLAPACK/blas_$PLAT.a $CLAPACK/F2CLIBS/libf2c.a

.o: %.c %.h
\t\$(CC) \$(CFLAGS) -c \$< -o \$@

clean:
\trm -f *.o *.so" > RPHAST/src/Makefile
else
cat <<EOF > RPHAST/src/Makevars
PKG_CFLAGS = -Iinclude/ -fgnu89-inline -DPHAST_VERSION=\"$PHAST_VERSION\" -DRPHAST -I$CLAPACK/INCLUDE
PKG_LIBS = \$(LAPACK_LIBS) \$(BLAS_LIBS) \$(FLIBS) -pthread
CSOURCES = `echo $CSOURCES`
OBJECTS = \$(CSOURCES:.c=.o)
all: \$(OJBECTS)
EOF
fi

R CMD INSTALL RPHAST/
