SUB = lib exoniphy phastCons phastOdds phastMotif phyloFit util 

CDIR = ${PWD}
INSTALLDIR = /cluster/bin/phast

all: 
	mkdir -p ../bin ../lib
	for dir in $(SUB) ; do cd ${CDIR}/$$dir && ${MAKE} ; done

install:
	cp ../bin/* ${INSTALLDIR}

doc:
	doxygen

clean:
	@for dir in $(SUB) ; do cd ${CDIR}/$$dir && ${MAKE} clean ; done
	rm -f ../bin/* ../lib/*
