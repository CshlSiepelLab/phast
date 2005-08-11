SUB = lib exoniphy phastCons phastOdds phastMotif phyloFit phyloBoot phyloP util 

CDIR = ${PWD}
INSTALLDIR = /cluster/bin/phast/${MACHTYPE}

all: 
	mkdir -p ../bin ../lib
	for dir in $(SUB) ; do cd ${CDIR}/$$dir && ${MAKE} ; done

install:
	@if [ "${MACHTYPE}" == "" ] ; then echo "ERROR: MACHTYPE must be defined." ; exit 1 ; fi
	mkdir -p ${INSTALLDIR}
	cp ../bin/* ${INSTALLDIR}

doc:
	doxygen

clean:
	@for dir in $(SUB) ; do cd ${CDIR}/$$dir && ${MAKE} clean ; done
	rm -f ../bin/* ../lib/*
